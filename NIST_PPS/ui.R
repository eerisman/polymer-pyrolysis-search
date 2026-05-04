library(shinyFiles)
library(plotly)
ui <- fluidPage(
  titlePanel(windowTitle = 'NIST Polymer Pyrolysis Search',
    fluidRow(
      column(8,h1('NIST Polymer Pyrolysis Search'), align = 'center'),
      column(4,img(height = 100, src='logo.png', align = 'center')))),
  
    fluidRow(
      column(2,
          shinyFilesButton('file', '*.d or netcdf file select', 'Please select a data file', FALSE),
          textOutput("path"),
          actionButton("RI", "Create RI calibration", class = "btn-success"),
          sliderInput('mfc', "Match Factor Cutoff", 700, min = 500, max = 999, step = 5 ),
          sliderInput('rip', "RI Penalty Rate", 50, min = 10, max = 1000, step = 5 ),
          sliderInput('snc', "S/N Cutoff", 50, min = 3, max = 1000),
          textOutput("system"),
          textOutput("system2"),
          actionButton("search", "Polymer Search", class = "btn-success"),
          textOutput("system3")
        ), 
      column(5,
          plotlyOutput("chromatogram.plotly")
        ),
      column(5,
          plotlyOutput("spectra.plotly")
                  )
  ),
  
  fluidRow(
    column(4, uiOutput("polymer.struct")),
    column(8, uiOutput("pyrolyzate.struct"))       
  ),
  
  fluidRow(
      column(4, DT::dataTableOutput("polymer")),
      column(8, DT::dataTableOutput("peaks"))
        )
  )