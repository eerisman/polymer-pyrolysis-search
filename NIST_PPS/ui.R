library(shinyFiles)

ui <- fluidPage(
  titlePanel("NIST Pyrolysis Polymer Search"),
  fluidRow(
      column(2,
          shinyFilesButton('file', '*.d or netcdf file select', 'Please select a data file', FALSE),
          textOutput("path"),
          actionButton("RI", "Create RI calibration", class = "btn-success"),
          tableOutput("table"),
          textOutput("system"),
          textOutput("system2"),
          actionButton("search", "Polymer Search", class = "btn-success"),
          textOutput("system3")
        ), 
      column(5,
          plotOutput("chromatogram")
        ),
      column(5,
          plotOutput("spectra")
        )
  ),
  
  fluidRow(
      column(4, DT::dataTableOutput("polymer")),
      column(8, DT::dataTableOutput("peaks"))
        )
)