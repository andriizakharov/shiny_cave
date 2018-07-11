#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("TEST"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
          sliderInput("control",
                      "This controls input",
                      min = 0, max = 10, value = 0),
         numericInput("num",
                     "Slope-slope correlation:",
                     min = 0,
                     max = 10,
                     value = 5)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("plot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    observeEvent(input$control, {
        updateNumericInput(session, "num", value = input$control)
    })
    
    out <- reactive(input$num)
    
    output$plot <- renderPlot({
        plot(0:10, rep(out(), 11), type = "l")
    })
}


# Run the application 
shinyApp(ui, server)

