library(bio3d)
library(shiny)
library(DT)

shinyUI(fluidPage(
  titlePanel("Kinesin Consensus Numbering Widget"),
  
  fluidRow(
    column(4,
           sliderInput("decimal", "Threshold:", 
                       min = 0, max = 1, value = 0.4, step= 0.1),
           tableOutput("values"))
  ),
  fluidRow(
    column(4,
           textInput("text", label = strong("Enter Kinesin ID")))
  ),
  wellPanel(
    div(style="display:inline-block",  
        h4("Results 1:"),
        dataTableOutput("output_table"))
  ),

  
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
  )
  
  ))

