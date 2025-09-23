source("./functions/elu2mspec.R", local = TRUE)
source("./functions/rical.R", local = TRUE)
source("./functions/mssearch.R", local = TRUE)

library(mssearchr)
library(mzR)
library(plotly)

server <- function(input, output, session) {
  volumes = getVolumes()
  shinyFileChoose(input, 'file', roots=volumes, filetypes=c('', 'd', 'cdf'))
  #create AMDIS onsite.ini file
  loc <- getwd()
  onsite <- readLines(paste0(loc,"/AMDIS/onsiteminimum.ini"))
  lib <- c(onsite, paste0("CALIBLIB=",loc,"/AMDIS/LIB/Local.cal"), paste0("LIB=",loc,"/AMDIS/LIB/ONSITE.MSL"))
  lib <- (gsub("/","\\\\",lib))
  writeLines(lib, "./AMDIS/onsite.ini")
  
  datapath <- reactive({parseFilePaths(volumes, input$file)$datapath})
    
    #display data path
    output$path <- renderText({datapath()})
    output$system <- renderText({if (length(datapath()) < 1){NULL
    }else if (grepl(".*\\.D", datapath(), ignore.case = T) & 
                !file.exists(paste0(datapath(),"/tic_front.csv"))){"Tic not found"}
          })
      
    #get tic data
    tic <- reactive({
            if (length(datapath()) < 1) {NULL
                }else if (file.exists(paste0(datapath(),"/tic_front.csv"))){
                read.csv(paste0(datapath(),"/tic_front.csv"), skip = 2, header = F, col.names = c("time", "abund"))
                }else if (grepl(".*\\.CDF", datapath(), ignore.case = T)){
                data.header <- header(openMSfile(datapath()))
                cbind(data.header$retentionTime/60, data.header$totIonCurrent)}
           })
    
    #checking for RI cal file
    rifile <- reactiveVal(
      if (file.exists("./AMDIS/LIB/Local.cal")){"RI cal file found"}
      else{"RI cal file not found"})
    
    #RI cal data for annotations modified
    ri.cal <- reactiveVal(if (file.exists("./AMDIS/LIB/Local.cal")){read.csv("./AMDIS/LIB/Local.cal", sep = "", header = F)}else{NULL})
    
    #run AMDIS for RI
    observeEvent(input$RI, {
      command <- if (grepl(".*\\.d", datapath(), ignore.case = T)){paste0("./AMDIS/AMDIS_32.exe ", datapath(),"\\data.ms")
        }else{paste0("./AMDIS/AMDIS_32.exe ", datapath())} #modified
      command <- paste0(gsub("/","\\\\",command), " \\/e")
      system(command)
      output$system3 <- renderText("AMDIS complete")
      if (file.exists(gsub("\\..*",".ELU", datapath()))){elupath <- gsub("\\..*",".ELU", datapath())}else{output$system2 <- renderText("No ELU")}
      elu2mspec(elupath)
      cal <- rical(paste0("./data/",gsub(".*/(.+).ELU", "\\1", elupath), ".MSPEC"))
      rifile(if (file.exists("./AMDIS/LIB/Local.cal")){"RI file found"}else{"not found"}) #updating the RI message
      ri.cal(read.csv("./AMDIS/LIB/Local.cal", sep = "", header = F)) #updating the ri calibration data for display on chromatogram
      })

    #plot tic
    output$chromatogram.plotly <-renderPlotly({
       if (is.null(tic())){plotly_empty(type = "scatter", mode = "lines")}else
       {layout(plot_ly(x=tic()[,1],y=tic()[,2], type = "scatter", mode = "lines"),
              xaxis = list(title = "time"),
              yaxis = list(title = "abundance"),
              title = gsub(".*/(.+)\\..*", "\\1", datapath()),
              annotations = lapply(1:nrow(ri.cal()), function(i){
                list(x = ri.cal()[i,1],
                     y = 1,
                     text = paste0("C",ri.cal()[i,2]/100),
                     showarrow = F,
                     xref = "x",
                     yref = "paper",
                     font = list(size = 15, weight = 'bold', color = 'red'))
              })
              )}

    })
    
    #display RI found message
    output$system2 <- renderText({rifile()})
    
    #run AMDIS on data
    observeEvent(input$search,{
      command1 <-if (grepl(".*\\.d", datapath(), ignore.case = T)){paste0("./AMDIS/AMDIS_32.exe ", datapath(),"\\data.ms")
      }else{paste0("./AMDIS/AMDIS_32.exe ", datapath())}#modified
      command1 <- paste0(gsub("/","\\\\",command1), " \\/x \\/e")
      system(command1)
      elupath <- if (file.exists(gsub("\\..*",".ELU", datapath()))){gsub("\\..*",".ELU", datapath())} #modified
      elu2mspec(elupath)
      results <- mssearch(paste0("./data/", gsub(".*/(.+)\\..*", "\\1", datapath()),".MSPEC"), input$snc, input$mfc, input$rip)  #modified
      output$polymer <- DT::renderDataTable(results[[1]], selection = 'single', rownames = FALSE)
      output$peaks <- DT::renderDataTable(results[[2]][,c(13,12,1,2,4:6,8:11,14)], selection = 'single', rownames = FALSE)
      
      output$spectra.plotly <- renderPlotly({
         s <- input$peaks_rows_selected
         if (is.null(s)){
           plotly_empty(type = "scatter", mode = "lines")
           return()}
         lib.entry <- results[[2]][s,]$idx
         query.entry <- results[[2]][s,]$qindx
           if (exists("lib.entry")&exists("query.entry")){
           plot_ly(type = "scatter", mode = "lines") %>%
               add_trace(x = as.vector(sapply(results[[3]][[query.entry]]$mz, function(x) c(x,x,NA))), 
                         y = as.vector(sapply(results[[3]][[query.entry]]$intst, function(y) c(0,y,NA))), 
                         name = "Query", line = list(color = "blue", width = 3)) %>%
               add_trace(x = as.vector(sapply(results[[4]][[lib.entry]]$mz, function(x) c(x,x,NA))),
                         y = -as.vector(sapply(results[[4]][[lib.entry]]$intst, function(y) c(0,y,NA))),
                         name = "Library", line = list(color = "red", width = 3)) %>%
           layout(
               title = "Spectra",
               xaxis = list(title = "m/z", 
                            range = c(min(c(results[[3]][[query.entry]]$mz,results[[4]][[lib.entry]]$mz))-5
                                      ,max(c(results[[3]][[query.entry]]$mz,results[[4]][[lib.entry]]$mz))+5)),
               yaxis = list(title = "rel. abund", range = c(-1000,1000))
           )   
             }
       }
       )
      })
    
}