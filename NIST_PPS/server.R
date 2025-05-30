source("./functions/elu2mspec.R", local = TRUE)
source("./functions/rical.R", local = TRUE)
source("./functions/mssearch.R", local = TRUE)

library(mssearchr)
library(mzR)

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
                read.csv(paste0(datapath(),"/tic_front.csv"), skip = 2)
                }else if (grepl(".*\\.CDF", datapath(), ignore.case = T)){
                data.header <- header(openMSfile(datapath()))
                cbind(data.header$retentionTime/60, data.header$totIonCurrent)}
           })
    
    #checking for RI cal file
    rifile <- reactiveVal(
      if (file.exists("./AMDIS/LIB/Local.cal")){"RI cal file found"}
      else{"RI cal file not found"})
    
    #run AMDIS for RI
    observeEvent(input$RI, {
      command <- if (grepl(".*\\.d", datapath(), ignore.case = T)){paste0("./AMDIS/AMDIS_32.exe ", datapath(),"\\data.ms")
        }else{paste0("./AMDIS/AMDIS_32.exe ", datapath())} #modified
      command <- paste0(gsub("/","\\\\",command), " \\/e")
      system(command)
      output$system2 <- renderText("AMDIS complete")
      if (file.exists(gsub("\\..*",".ELU", datapath()))){elupath <- gsub("\\..*",".ELU", datapath())}else{output$system2 <- renderText("No ELU")}
      elu2mspec(elupath)
      cal <- rical(paste0("./data/",gsub(".*/(.+).ELU", "\\1", elupath), ".MSPEC"))
      output$chromatogram <-renderPlot(
        if (is.null(tic())){plot.new()}else
        {plot(tic(), xlab = "time", ylab = "abundance", type = "l", main = gsub(".*/(.+)\\..*", "\\1", datapath()))  #modified
        text(cal[,1], max(tic()[,2]), paste0("C",cal[,2]/100), col = "RED")
     })
    })

    #updating that RI cal file found
    observeEvent(input$RI, {
      rifile(if (file.exists("./AMDIS/LIB/Local.cal")){"RI file found"}else{"not found"})
    })
    
    #plot tic  
    output$chromatogram <-renderPlot({
      if (is.null(tic())){plot.new()}else{
        plot(tic(), xlab = "time", ylab = "abundance", type = "l", main = gsub(".*/(.+)\\..*", "\\1", datapath()))  #modified 
        }
      })
    
    #display RI found message
    output$system2 <- renderText({rifile()})
    
    #run AMDIS on data ##needs work
    observeEvent(input$search,{
      command1 <-if (grepl(".*\\.d", datapath(), ignore.case = T)){paste0("./AMDIS/AMDIS_32.exe ", datapath(),"\\data.ms")
      }else{paste0("./AMDIS/AMDIS_32.exe ", datapath())}#modified
      command1 <- paste0(gsub("/","\\\\",command1), " \\/x \\/e")
      system(command1)
      elupath <- if (file.exists(gsub("\\..*",".ELU", datapath()))){gsub("\\..*",".ELU", datapath())} #modified
      elu2mspec(elupath)
      results <- mssearch(paste0("./data/", gsub(".*/(.+)\\..*", "\\1", datapath()),".MSPEC"))  #modified
      output$polymer <- DT::renderDataTable(results[[1]], rownames = FALSE)
      output$peaks <- DT::renderDataTable(results[[2]][,c(13,12,1,2,4:6,8:11,14)], selection = 'single', rownames = FALSE)
      
      output$spectra <- renderPlot({
        s <- input$peaks_rows_selected
        if (is.null(s)){
          plot.new()
          return()}
        lib.entry <- results[[2]][s,]$idx
        query.entry <- results[[2]][s,]$qindx
          if (exists("lib.entry")&exists("query.entry")){
          plot(results[[4]][[lib.entry]]$mz,
                 results[[4]][[lib.entry]]$intst,
                 xlab = "m/z", ylab = "abundance", ylim = c(-1000,1000), 
                 type = "h", main = "Spectra", lwd = 2, col = "blue")
          lines(results[[3]][[query.entry]]$mz,
                -results[[3]][[query.entry]]$intst,
                type = "h", lwd = 2, col = "red")
          mtext("Query", side = 3, col = "blue")
          mtext("Library", side = 1, col = "red")
            }
      }
      )
      })
    
}