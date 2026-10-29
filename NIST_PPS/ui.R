library(shinyFiles)
library(plotly)

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .well { padding: 12px; margin-bottom: 12px; background-color: #f8f9fa; border: 1px solid #e9ecef; border-radius: 6px; }
      /* Ensure DataTables container and tables stay within card boundaries and support horizontal scroll */
      .dataTables_wrapper { font-size: 13px; overflow-x: auto; max-width: 100%; }
      table.dataTable { width: 100% !important; }
      h4, h5 { color: #2c3e50; font-weight: 600; margin-top: 2px; }
      .btn-custom { font-weight: bold; width: 100%; }
      /* Ensure SVG structure images resize cleanly and never extend off screen */
      .struct-container img { max-width: 100%; height: auto; max-height: 120px; display: block; margin: 0 auto; }
      .struct-box { min-height: 110px; display: flex; align-items: center; justify-content: center; }
      .nav-tabs { margin-bottom: 15px; font-weight: bold; font-size: 15px; }
      .help-section { line-height: 1.6; font-size: 14px; color: #2c3e50; }
      .help-card { background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 18px; margin-bottom: 16px; box-shadow: 0 1px 3px rgba(0,0,0,0.05); }
      .help-card h4 { color: #1a202c; font-weight: 700; margin-top: 0; margin-bottom: 12px; border-bottom: 2px solid #edf2f7; padding-bottom: 6px; }
      .help-card h5 { color: #2b6cb0; font-weight: 600; margin-top: 10px; margin-bottom: 6px; }
      .badge-param { background-color: #edf2f7; color: #2d3748; padding: 2px 6px; border-radius: 4px; font-family: monospace; font-weight: bold; }
    "))
  ),

  # Header Bar
  fluidRow(style = "background-color: #2c3e50; color: white; padding: 10px 20px; margin-bottom: 15px; border-radius: 4px;",
    column(8,
      div(style = "display: flex; align-items: center; gap: 12px;",
        h2("NIST Polymer Pyrolysis Search", style = "margin: 0; color: white; display: inline-block;"),
        span("v1.0.0", style = "font-size: 14px; background-color: #34495e; padding: 3px 8px; border-radius: 4px; color: #ecf0f1; font-weight: bold; font-family: monospace;")
      )
    ),
    column(4, img(height = 45, src = "logo.svg", align = "right"))
  ),

  # Main Content Tabs
  tabsetPanel(id = "main_tabs",
    tabPanel("Dashboard",
      br(),
      fluidRow(
        # --- LEFT COLUMN: CONTROL PANEL ---
        column(2,
          div(class = "well",
            h4("Data File"),
            shinyFilesButton('file', 'Select File (*.d / *.cdf)', 'Please select a data file', FALSE, class = "btn-primary btn-sm btn-custom"),
            div(style = "font-size: 11px; color: #555; word-break: break-all; margin-top: 5px;", textOutput("path")),
            textOutput("tic_status")
          ),
          div(class = "well",
            h4("RI Calibration"),
            actionButton("RI", "Create RI Calibration", class = "btn-info btn-sm btn-custom"),
            div(style = "font-size: 11px; font-weight: bold; margin-top: 4px;", textOutput("ri_status"))
          ),
          div(class = "well",
            h4("AMDIS Deconvolution"),
            actionButton("amdis", "Run AMDIS", class = "btn-warning btn-sm btn-custom"),
            div(style = "font-size: 11px; font-weight: bold; margin-top: 4px;", textOutput("amdis_status"))
          ),
          div(class = "well",
            h4("Parameters"),
            sliderInput('snc', "S/N Cutoff", 50, min = 3, max = 1000, ticks = FALSE),
            sliderInput('mfc', "Match Factor Cutoff", 700, min = 500, max = 999, step = 5, ticks = FALSE),
            sliderInput('rip', "RI Penalty Rate", 50, min = 10, max = 1000, step = 5, ticks = FALSE),
            checkboxInput('waive_out_of_range', "Waive RI penalty outside calibration range", value = TRUE)
          ),
          div(class = "well",
            actionButton("search", "Polymer Search", class = "btn-success btn-custom"),
            div(style = "font-size: 11px; font-weight: bold; margin-top: 4px;", textOutput("search_status"))
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
    ),

    tabPanel("Help & Instructions",
      br(),
      div(class = "help-section",
        fluidRow(
          column(12,
            div(class = "help-card",
              h3("NIST Polymer Pyrolysis Search (NPPS) v1.0.0", style = "margin-top: 0; color: #1a202c;"),
              p("The NIST Polymer Pyrolysis Search application is designed for analyzing Pyrolysis Gas Chromatography / Mass Spectrometry (Py-GC/MS) data to search, identify, and quantify polymers and pyrolyzates in query samples.")
            )
          )
        ),
        fluidRow(
          column(6,
            div(class = "help-card",
              h4("Workflow & Step-by-Step Instructions"),
              ol(style = "padding-left: 20px;",
                tags$li(
                  tags$strong("Select Data File: "),
                  "Click ", tags$span(class = "badge-param", "Select File"), " to choose a Py-GC/MS data folder/file (Agilent ", tags$code(".d"), " or NetCDF ", tags$code(".cdf"), "). Note: Data file paths should not contain spaces."
                ),
                tags$li(
                  tags$strong("Create RI Calibration: "),
                  "Select alkane standard data and click ", tags$span(class = "badge-param", "Create RI Calibration"), ". This runs AMDIS to build/update the local Retention Index calibration file (", tags$code("./AMDIS/LIB/Local.cal"), ")."
                ),
                tags$li(
                  tags$strong("Run AMDIS Deconvolution: "),
                  "Click ", tags$span(class = "badge-param", "Run AMDIS"), " to deconvolute peaks from the query pyrogram and extract individual spectra into a combined ", tags$code(".MSPEC"), " file. Polymer Search automatically executes upon completion."
                ),
                tags$li(
                  tags$strong("Adjust Parameters & Re-Search: "),
                  "Modify search parameter sliders or options as needed and click ", tags$span(class = "badge-param", "Polymer Search"), ". Re-searching uses cached ", tags$code(".MSPEC"), " data without requiring AMDIS to re-run."
                )
              )
            ),
            div(class = "help-card",
              h4("Interactive Dashboard Features"),
              ul(style = "padding-left: 20px;",
                tags$li(
                  tags$strong("Chromatogram (TIC) & RI Markers: "),
                  "Displays Total Ion Chromatogram with Retention Index calibration marker overlay lines. Selected pyrolyzate peaks display an arrow marker at their retention time."
                ),
                tags$li(
                  tags$strong("Mass Spectral Head-to-Tail Comparison: "),
                  "Selecting a pyrolyzate peak row dynamically plots the query spectrum (top/blue) against the reference library spectrum (bottom/red)."
                ),
                tags$li(
                  tags$strong("Polymer & Pyrolyzate Chemical Structures: "),
                  "Displays chemical structure diagrams for the selected polymer hit and pyrolyzate peak."
                ),
                tags$li(
                  tags$strong("Expected Pyrolyzates Table: "),
                  "Selecting a polymer candidate updates the expected pyrolyzates table showing theoretical marker compounds for that polymer."
                )
              )
            )
          ),
          column(6,
            div(class = "help-card",
              h4("Search Parameter Reference"),
              ul(style = "padding-left: 20px;",
                tags$li(
                  tags$strong("S/N Cutoff (snc): "),
                  "Filters out components with Signal-to-Noise ratio below the selected cutoff (Range: 3 – 1000, Default: 50)."
                ),
                tags$li(
                  tags$strong("Match Factor Cutoff (mfc): "),
                  "Minimum required mass spectral similarity score (Range: 500 – 999, Default: 700)."
                ),
                tags$li(
                  tags$strong("RI Penalty Rate (rip): "),
                  "Controls penalty applied to the match factor when query RI deviates from library RI (Range: 10 – 1000, Default: 50)."
                ),
                tags$li(
                  tags$strong("Waive RI penalty outside calibration range: "),
                  "When checked (default), peaks with retention index values outside the calibrated range in ", tags$code("Local.cal"), " are exempted from RI deviation penalties."
                )
              )
            ),
            div(class = "help-card",
              h4("Exported Search Results & File Location"),
              p("Each search automatically creates a unique parameter-and-timestamp-named subfolder inside ", tags$code("./results/"), " containing:"),
              ul(style = "padding-left: 20px;",
                tags$li(tags$code("*_search_parameters_*.txt"), ": Summary of search parameters, date/time, and application version (", tags$code("v1.0.0"), ")."),
                tags$li(tags$code("*_rawsearchresults_*.csv"), ": Unfiltered raw pyrolyzate library hits."),
                tags$li(tags$code("*_filteredsearchresults_*.csv"), ": Filtered hitlist sorted by modified match factor."),
                tags$li(tags$code("*_polymerlist_*.csv"), ": Final polymer identification results sorted by final score.")
              )
            )
          )
        )
      )
    )
  )
)
