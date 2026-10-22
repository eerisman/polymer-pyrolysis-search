source("./functions/elu2mspec.R", local = TRUE)
source("./functions/rical.R", local = TRUE)
source("./functions/mssearch.R", local = TRUE)

library(mssearchr)
library(mzR)
library(plotly)

server <- function(input, output, session) {
  volumes <- getVolumes()()
  shinyFileChoose(input, 'file', roots = volumes, filetypes = c('', 'd', 'cdf'))
  
  # create AMDIS onsite.ini file
  loc <- getwd()
  if (file.exists(paste0(loc, "/AMDIS/onsiteminimum.ini"))) {
    onsite <- readLines(paste0(loc, "/AMDIS/onsiteminimum.ini"))
    lib <- c(onsite, paste0("CALIBLIB=", loc, "/AMDIS/LIB/Local.cal"), paste0("LIB=", loc, "/AMDIS/LIB/ONSITE.MSL"))
    lib <- gsub("/", "\\\\", lib)
    writeLines(lib, "./AMDIS/onsite.ini")
  }
  
  # results file reactive
  results <- reactiveVal(NULL)
  amdis_status <- reactiveVal("AMDIS not run")
  output$amdis_status <- renderText({ amdis_status() })
  
  datapath <- reactive({
    results(NULL)
    dp <- parseFilePaths(volumes, input$file)$datapath
    if (length(dp) > 0 && !is.na(dp) && dp != "") {
      clean_name <- gsub(" ", "_", tools::file_path_sans_ext(basename(dp)))
      mspec_path <- paste0("./data/", clean_name, ".MSPEC")
      if (file.exists(mspec_path)) {
        amdis_status("AMDIS data ready")
      } else {
        amdis_status("AMDIS not run")
      }
    } else {
      amdis_status("No file selected")
    }
    dp
  })
    
  # display data path
  output$path <- renderText({ datapath() })
  output$system <- renderText({
    if (length(datapath()) < 1 || is.na(datapath()) || datapath() == "") {
      NULL
    } else if (grepl(".*\\.D", datapath(), ignore.case = TRUE) &
               !file.exists(paste0(datapath(), "/tic_front.csv")) &
               !file.exists(paste0(datapath(), "/tic_back.csv"))) {
      "Tic not found"
    }
  })
      
  # get tic data
  tic <- reactive({
    if (length(datapath()) < 1 || is.na(datapath()) || datapath() == "") {
      NULL
    } else if (file.exists(paste0(datapath(), "/tic_front.csv"))) {
      read.csv(paste0(datapath(), "/tic_front.csv"), skip = 2, header = FALSE, col.names = c("time", "abund"))
    } else if (file.exists(paste0(datapath(), "/tic_back.csv"))) {
      read.csv(paste0(datapath(), "/tic_back.csv"), skip = 2, header = FALSE, col.names = c("time", "abund"))
    } else if (grepl(".*\\.CDF", datapath(), ignore.case = TRUE)) {
      msfile <- openMSfile(datapath())
      data.header <- header(msfile)
      close(msfile)
      cbind(data.header$retentionTime / 60, data.header$totIonCurrent)
    }
  })
    
  # checking for RI cal file
  rifile <- reactiveVal(
    if (file.exists("./AMDIS/LIB/Local.cal")) "RI cal file found" else "RI cal file not found"
  )
    
  # RI cal data for annotations modified
  ri.cal <- reactiveVal(
    if (file.exists("./AMDIS/LIB/Local.cal")) read.csv("./AMDIS/LIB/Local.cal", sep = "", header = FALSE) else NULL
  )
    
  # run AMDIS for RI
  observeEvent(input$RI, {
    req(datapath())
    if (length(datapath()) < 1 || is.na(datapath()) || datapath() == "") return()

    command <- if (grepl(".*\\.d", datapath(), ignore.case = TRUE)) {
      paste0("./AMDIS/AMDIS_32.exe ", datapath(), "\\data.ms")
    } else {
      paste0("./AMDIS/AMDIS_32.exe ", datapath())
    }
    command <- paste0(gsub("/", "\\\\", command), " \\/e")
    system(command)
    output$system3 <- renderText("AMDIS complete")

    elupath <- sub("\\.[^/\\\\]+$", ".ELU", datapath())
    if (file.exists(elupath)) {
      elu2mspec(elupath)
      clean_name <- gsub(" ", "_", tools::file_path_sans_ext(basename(elupath)))
      cal <- rical(paste0("./data/", clean_name, ".MSPEC"))
      rifile(if (file.exists("./AMDIS/LIB/Local.cal")) "RI file found" else "not found")
      if (file.exists("./AMDIS/LIB/Local.cal")) {
        ri.cal(read.csv("./AMDIS/LIB/Local.cal", sep = "", header = FALSE))
      }
    } else {
      output$system2 <- renderText("No ELU")
    }
  })

  # plot tic
  output$chromatogram.plotly <- renderPlotly({
    if (is.null(tic())) {
      return(plotly_empty(type = "scatter", mode = "lines"))
    }
    
    chrom <- plot_ly(x = tic()[,1], y = tic()[,2], type = "scatter", mode = "lines") %>%
      layout(
        xaxis = list(title = "Time (min)"),
        yaxis = list(title = "Abundance"),
        title = tools::file_path_sans_ext(basename(datapath()))
      )
    
    ann.ri <- NULL
    shapes.ri <- NULL
    if (!is.null(ri.cal()) && nrow(ri.cal()) > 0) {
      y_tiers <- c(1.02, 0.91, 0.80)
      ann.ri <- lapply(seq_len(nrow(ri.cal())), function(i) {
        y_pos <- y_tiers[((i - 1) %% length(y_tiers)) + 1]
        list(
          x = ri.cal()[i, 1],
          y = y_pos,
          text = paste0("C", ri.cal()[i, 2] / 100),
          showarrow = FALSE,
          xref = "x",
          yref = "paper",
          font = list(size = 12, weight = 'bold', color = 'red')
        )
      })

      shapes.ri <- lapply(seq_len(nrow(ri.cal())), function(i) {
        list(
          type = "line",
          x0 = ri.cal()[i, 1],
          x1 = ri.cal()[i, 1],
          y0 = 0,
          y1 = 1,
          xref = "x",
          yref = "paper",
          line = list(color = "rgba(255, 0, 0, 0.35)", width = 1, dash = "dot")
        )
      })
    }

    ann.marker <- NULL
    if (!is.null(results()) && !is.null(input$peaks_rows_selected) && length(input$peaks_rows_selected) > 0) {
      peak_row <- input$peaks_rows_selected
      if (peak_row <= nrow(results()[[2]])) {
        rt_val <- results()[[2]][peak_row, 12]
        if (!is.na(rt_val)) {
          closest_idx <- which.min(abs(tic()[, 1] - rt_val))
          ann.marker <- list(list(
            x = rt_val,
            y = tic()[closest_idx, 2],
            text = "",
            ax = 0,
            ay = -30,
            arrowcolor = 'green',
            arrowwidth = 3
          ))
        }
      }
    }

    chrom %>% layout(annotations = c(ann.ri, ann.marker), shapes = shapes.ri)
  })

  # display RI found message
  output$system2 <- renderText({ rifile() })

  # Function to run AMDIS deconvolution
  run_amdis_deconv <- function() {
    req(datapath())
    if (length(datapath()) < 1 || is.na(datapath()) || datapath() == "") return(FALSE)

    command1 <- if (grepl(".*\\.d", datapath(), ignore.case = TRUE)) {
      paste0("./AMDIS/AMDIS_32.exe ", datapath(), "\\data.ms")
    } else {
      paste0("./AMDIS/AMDIS_32.exe ", datapath())
    }
    command1 <- paste0(gsub("/", "\\\\", command1), " \\/x \\/e")
    system(command1)

    elupath <- sub("\\.[^/\\\\]+$", ".ELU", datapath())
    if (file.exists(elupath)) {
      elu2mspec(elupath)
      amdis_status("AMDIS complete")
      return(TRUE)
    } else {
      amdis_status("No ELU file generated")
      return(FALSE)
    }
  }

  # Dedicated AMDIS button action
  observeEvent(input$amdis, {
    run_amdis_deconv()
  })

  # Polymer search on data (uses cached MSPEC if available)
  observeEvent(input$search, {
    req(datapath())
    if (length(datapath()) < 1 || is.na(datapath()) || datapath() == "") return()

    clean_name <- gsub(" ", "_", tools::file_path_sans_ext(basename(datapath())))
    mspec_path <- paste0("./data/", clean_name, ".MSPEC")

    if (!file.exists(mspec_path)) {
      amdis_ok <- run_amdis_deconv()
      if (!amdis_ok || !file.exists(mspec_path)) {
        output$system3 <- renderText("AMDIS deconvolution failed")
        return()
      }
    }

    res <- mssearch(mspec_path, input$snc, input$mfc, input$rip)
    results(res)
    output$system3 <- renderText("Search complete")
  })

  # Top-level reactive UI and outputs
  output$polymer <- DT::renderDataTable({
    req(results())
    results()[[1]]
  }, selection = 'single', rownames = FALSE, options = list(scrollX = TRUE, autoWidth = TRUE))

  output$peaks <- DT::renderDataTable({
    req(results())
    df <- results()[[2]]
    if (nrow(df) > 0 && ncol(df) >= 13) {
      df[, c(12, 11, 1, 2, 4:5, 7:10, 13)]
    } else {
      df
    }
  }, selection = 'single', rownames = FALSE, options = list(scrollX = TRUE, autoWidth = TRUE))

  output$Search.polymer <- renderUI({
    req(results())
    tags$h2("Polymer Search Results", style = "text-align: center;")
  })

  output$Search.pyrolyzate <- renderUI({
    req(results())
    tags$h2("Pyrolyzate Search Results", style = "text-align: center;")
  })

  output$polymer.struct <- renderUI({
    req(results(), input$polymer_rows_selected)
    p_row <- input$polymer_rows_selected
    polym_df <- results()[[1]]
    if (p_row <= nrow(polym_df)) {
      polymer.name <- sub(" \\(.*\\)", "", polym_df[p_row, 1])
      polymer.path <- paste0("polymers/", polymer.name, ".svg")
      if (file.exists(paste0("www/", polymer.path))) {
        img(src = polymer.path, height = 100)
      } else {
        div(style = "color: #888; font-style: italic;", "Structure not available")
      }
    }
  })

  output$pyrolyzate.struct <- renderUI({
    req(results(), input$peaks_rows_selected)
    pk_row <- input$peaks_rows_selected
    peaks_df <- results()[[2]]
    if (pk_row <= nrow(peaks_df)) {
      pyrolyzate.index <- peaks_df[pk_row, 3]
      pyrolyzate.path <- paste0("pyrolyzates/", pyrolyzate.index, ".svg")
      if (file.exists(paste0("www/", pyrolyzate.path))) {
        img(src = pyrolyzate.path, height = 100)
      } else {
        div(style = "color: #888; font-style: italic;", "Structure not available")
      }
    }
  })

  output$expected <- renderUI({
    req(results(), input$polymer_rows_selected)
    h2("Expected Pyrolyzates", style = "text-align: center;")
  })

  output$pyrolyzates.expected <- DT::renderDT({
    req(results(), input$polymer_rows_selected)
    p_row <- input$polymer_rows_selected
    polym_df <- results()[[1]]
    if (p_row > nrow(polym_df)) return(NULL)

    polymer.name <- sub(" \\(.*\\)", "", polym_df[p_row, 1])
    lib_list <- results()[[4]]

    markers <- which(sapply(lib_list, function(l) {
      if (!is.null(l$polymer)) grepl(polymer.name, l$polymer, fixed = TRUE) else FALSE
    }))

    if (length(markers) > 0) {
      expected_df <- do.call(rbind, lapply(markers, function(m) {
        item <- lib_list[[m]]
        data.frame(
          Name = ifelse(is.null(item$name), "", item$name),
          Formula = ifelse(is.null(item$formula), "", item$formula),
          MW = ifelse(is.null(item$mw), "", item$mw),
          RI = ifelse(is.null(item$ri), "", item$ri),
          stringsAsFactors = FALSE
        )
      }))
      expected_df[order(expected_df$RI),]
    } else {
      NULL
    }
  }, selection = 'single', rownames = FALSE, options = list(scrollX = TRUE, autoWidth = TRUE))

  output$spectra.plotly <- renderPlotly({
    s <- input$peaks_rows_selected
    if (is.null(s) || is.null(results())) {
      return(plotly_empty(type = "scatter", mode = "lines"))
    }
    peaks_df <- results()[[2]]
    if (s > nrow(peaks_df)) {
      return(plotly_empty(type = "scatter", mode = "lines"))
    }

    lib.entry <- peaks_df[s, ]$idx
    query.entry <- peaks_df[s, ]$qindx

    query_list <- results()[[3]]
    lib_list <- results()[[4]]

    if (!is.null(query.entry) && !is.na(query.entry) && query.entry <= length(query_list) &&
        !is.null(lib.entry) && !is.na(lib.entry) && lib.entry <= length(lib_list)) {
      
      q_spec <- query_list[[query.entry]]
      l_spec <- lib_list[[lib.entry]]
      
      if (!is.null(q_spec$mz) && !is.null(l_spec$mz)) {
        q_x <- as.vector(sapply(q_spec$mz, function(x) c(x, x, NA)))
        q_y <- as.vector(sapply(q_spec$intst, function(y) c(0, y, NA)))
        l_x <- as.vector(sapply(l_spec$mz, function(x) c(x, x, NA)))
        l_y <- as.vector(sapply(l_spec$intst, function(y) c(0, -y, NA)))

        all_mz <- c(q_spec$mz, l_spec$mz)
        min_mz <- min(all_mz, na.rm = TRUE) - 5
        max_mz <- max(all_mz, na.rm = TRUE) + 5

        plot_ly(type = "scatter", mode = "lines") %>%
          add_trace(x = q_x, y = q_y, name = "Query", line = list(color = "blue", width = 3)) %>%
          add_trace(x = l_x, y = l_y, name = "Library", line = list(color = "red", width = 3)) %>%
          layout(
            title = "Spectra",
            xaxis = list(title = "m/z", range = c(min_mz, max_mz)),
            yaxis = list(title = "rel. abund", range = c(-1000, 1000))
          )
      } else {
        plotly_empty(type = "scatter", mode = "lines")
      }
    } else {
      plotly_empty(type = "scatter", mode = "lines")
    }
  })

  session$onSessionEnded(function() { stopApp() })
}
