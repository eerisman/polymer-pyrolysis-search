library(shinyFiles)
library(plotly)

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .well { padding: 12px; margin-bottom: 12px; background-color: #f8f9fa; border: 1px solid #e9ecef; border-radius: 6px; }
      .dataTables_wrapper { font-size: 13px; }
      h4, h5 { color: #2c3e50; font-weight: 600; margin-top: 2px; }
      .btn-custom { font-weight: bold; width: 100%; }
      /* Ensure SVG structure images resize cleanly and never extend off screen */
      .struct-container img { max-width: 100%; height: auto; max-height: 120px; display: block; margin: 0 auto; }
      .struct-box { min-height: 110px; display: flex; align-items: center; justify-content: center; }
    "))
  ),

  # Header Bar
  fluidRow(style = "background-color: #2c3e50; color: white; padding: 10px 20px; margin-bottom: 15px; border-radius: 4px;",
    column(8, h2("NIST Polymer Pyrolysis Search", style = "margin: 0; color: white;")),
    column(4, img(height = 45, src = "logo.png", align = "right"))
  ),

  fluidRow(
    # --- LEFT COLUMN: CONTROL PANEL ---
    column(2,
      div(class = "well",
        h4("Data File"),
        shinyFilesButton('file', 'Select File (*.d / *.cdf)', 'Please select a data file', FALSE, class = "btn-primary btn-sm btn-custom"),
        div(style = "font-size: 11px; color: #555; word-break: break-all; margin-top: 5px;", textOutput("path")),
        textOutput("system")
      ),
      div(class = "well",
        h4("RI Calibration"),
        actionButton("RI", "Create RI Calibration", class = "btn-info btn-sm btn-custom"),
        div(style = "font-size: 11px; font-weight: bold; margin-top: 4px;", textOutput("system2"))
      ),
      div(class = "well",
        h4("Parameters"),
        sliderInput('snc', "S/N Cutoff", 50, min = 3, max = 1000, ticks = FALSE),
        sliderInput('mfc', "Match Factor Cutoff", 700, min = 500, max = 999, step = 5, ticks = FALSE),
        sliderInput('rip', "RI Penalty Rate", 50, min = 10, max = 1000, step = 5, ticks = FALSE)
      ),
      div(class = "well",
        actionButton("search", "Polymer Search", class = "btn-success btn-custom"),
        textOutput("system3")
      )
    ),

    # --- RIGHT COLUMN: DASHBOARD CONTENT ---
    column(10,
      # ROW 1: INTERACTIVE PLOTS
      fluidRow(
        column(6,
          div(class = "well",
            h5("Chromatogram (TIC)"),
            plotlyOutput("chromatogram.plotly", height = "260px")
          )
        ),
        column(6,
          div(class = "well",
            h5("Mass Spectra Comparison (Query vs Library)"),
            plotlyOutput("spectra.plotly", height = "260px")
          )
        )
      ),

      # ROW 2: DEDICATED STRUCTURES ROW
      fluidRow(
        column(4,
          div(class = "well struct-box struct-container",
            div(style = "text-align: center; width: 100%;",
              h5("Polymer Structure", style = "margin-bottom: 8px;"),
              uiOutput("polymer.struct")
            )
          )
        ),
        column(8,
          div(class = "well struct-box struct-container",
            div(style = "text-align: center; width: 100%;",
              h5("Pyrolyzate Structure", style = "margin-bottom: 8px;"),
              uiOutput("pyrolyzate.struct")
            )
          )
        )
      ),

      # ROW 3: SEARCH RESULT TABLES
      fluidRow(
        column(4,
          div(class = "well",
            uiOutput("Search.polymer"),
            DT::dataTableOutput("polymer"),
            br(),
            uiOutput("expected"),
            DT::dataTableOutput("pyrolyzates.expected")
          )
        ),
        column(8,
          div(class = "well",
            uiOutput("Search.pyrolyzate"),
            DT::dataTableOutput("peaks")
          )
        )
      )
    )
  )
)
