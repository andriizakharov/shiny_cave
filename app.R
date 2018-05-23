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
   titlePanel("CAVE-SOS sim"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("cor",
                     "Slope-slope correlation:",
                     min = 0,
                     max = 1,
                     value = 0.5)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         htmlOutput("caption1"),
         htmlOutput("distText")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   output$caption1 <- renderText({"SOS:"})
   output$distText <- renderText({
      # generate bins based on input$bins from ui.R
        mean_icept_x <- 0
        mean_slope_x <- -14
        var_icept_x <- 100
        var_slope_x <- 200
        error_var_x <- 25
        cov_icept_x_slope_x <- 0
        
        mean_icept_y <- 0
        mean_slope_y <- -5
        var_icept_y <- 100
        var_slope_y <- 200
        error_var_y <- 25
        cov_icept_y_slope_y <- 0
        
        cov_icept_x_icept_y <- 60
        cov_icept_x_slope_y <- 0
        cov_icept_y_slope_x <- 0
        
        cor_slope_x_slope_y <- input$cor
        cov_slope_x_slope_y <- cor_slope_x_slope_y * sqrt(var_slope_x*var_slope_y)
        
        var_x <- var_icept_x + 1/3*var_slope_x + cov_icept_x_slope_x + 
           1/12*mean_slope_x**2 + error_var_x
        var_y <- var_icept_y + 1/3*var_slope_y + cov_icept_y_slope_y + 
           1/12*mean_slope_y**2 + error_var_y
        
        cov_x_y <- cov_icept_x_icept_y + 
           1/2*(cov_icept_x_slope_y + cov_icept_y_slope_x) +
           1/3*cov_slope_x_slope_y +
           1/12*mean_slope_x*mean_slope_y
        
        sos <- 1 - 12*var_y*(mean_slope_y*var_x - mean_slope_x*cov_x_y)**2 /
           ((var_y*var_x - cov_x_y**2) * (12*var_x - mean_slope_x**2)*mean_slope_y**2)

      # draw the histogram with the specified number of bins
        sos
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

