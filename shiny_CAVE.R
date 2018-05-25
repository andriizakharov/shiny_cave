library(semper)
source("bivariateLgcm.R")
source("toOpenMx.R")
source("input_presets.R")

library(shiny)
library(ggplot2)
library(tidyverse)

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("Bivariate LGCM simulation"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            titlePanel("Presets:"),
            
            selectInput("preset",
                        "Choose a preset:",
                        choices = list(
                            Preset_1 = preset1,
                            Preset_2 = preset2,
                            Preset_3 = preset3
                        )),
            
            titlePanel("Variable X:"),
            
            # Input params for variable X
            numericInput("mean_icept_x",
                         "X: intercept mean",
                         value = 0),
            numericInput("mean_slope_x",
                         "X: slope mean",
                         value = 0),
            numericInput("var_icept_x",
                         "X: intercept variance",
                         value = 0),
            numericInput("var_slope_x",
                         "X: slope variance",
                         value = 0),
            numericInput("error_var_x",
                         "X: error variance",
                         value = 0),
            numericInput("cov_icept_x_slope_x",
                         "X: intercept-slope covariance",
                         value = 0),
            
            titlePanel("Variable Y:"),
            
            # Input params for variable Y
            numericInput("mean_icept_y",
                         "Y: intercept mean",
                         value = 0),
            numericInput("mean_slope_y",
                         "Y: slope mean",
                         value = 0),
            numericInput("var_icept_y",
                         "Y: intercept variance",
                         value = 0),
            numericInput("var_slope_y",
                         "Y: slope variance",
                         value = 0),
            numericInput("error_var_y",
                         "Y: error variance",
                         value = 0),
            numericInput("cov_icept_y_slope_y",
                         "Y: intercept-slope covariance",
                         value = 0),
            
            titlePanel("Cross-variable:"),
            
            # Input intervariable params
            numericInput("cov_icept_x_icept_y",
                         "XY: intercept covariance",
                         value = 0),
            numericInput("cov_icept_x_slope_y",
                         "XY: covariance intercept X - slope Y",
                         value = 0),
            numericInput("cov_icept_y_slope_x",
                         "XY: covariance intercept Y - slope X",
                         value = 0),
            
            sliderInput("cor_slope_x_slope_y",
                        "XY: slope-slope correlation",
                        min = 0,
                        max = 1,
                        value = 0.5),
            
            titlePanel("Other:"),
            
            # How many timepoints should the models simulate?
            numericInput("timepoints",
                         "Number of timepoints per model",
                         min = 0,
                         max = 100,
                         value = 10),
            
            # How many simulated trajectories?
            numericInput("sim_sample",
                         "Number of simulated trajectories",
                         min = 0,
                         max = 100,
                         value = 20),
            
            # Button to simulate data
            actionButton("sim_data",
                         "Simulate!")
        ),
        
        # Show 
        mainPanel(
            titlePanel("SOS:"),
            htmlOutput("sos_val"),
            
            plotOutput("lgcm_graphs")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$sos_val <- renderText({
        
        # Define the models given all input parameters
        lgcm_x <- lgcm(
            timepoints = 1:input$timepoints,
            intercept.variance = input$var_icept_x,
            slope.variance = input$var_slope_x,
            residual.variance = input$error_var_x,
            intercept.slope.covariance = input$cov_icept_x_slope_x,
            intercept.mean = input$mean_icept_x,
            slope.mean = input$mean_slope_x
        )
        lgcm_y <- lgcm(
            timepoints = 1:input$timepoints,
            intercept.variance = input$var_icept_y,
            slope.variance = input$var_slope_y,
            residual.variance = input$error_var_y,
            intercept.slope.covariance = input$cov_icept_y_slope_y,
            intercept.mean = input$mean_icept_y,
            slope.mean = input$mean_slope_y
        )
        lgcm_xy <- bivariateLgcm(
            lgcm_x, lgcm_y, 
            icept.x.icept.y.covariance = input$cov_icept_x_icept_y, 
            slope.x.slope.y.covariance = input$cov_slope_x_slope_y
        )
        
        # Calculate some values
        cov_slope_x_slope_y <- input$cor_slope_x_slope_y * 
            sqrt(input$var_slope_x * input$var_slope_y)
        
        var_x <- input$var_icept_x + 1/3*input$var_slope_x + input$cov_icept_x_slope_x + 
            1/12*input$mean_slope_x**2 + input$error_var_x
        var_y <- input$var_icept_y + 1/3*input$var_slope_y + input$cov_icept_y_slope_y + 
            1/12*input$mean_slope_y**2 + input$error_var_y
        
        cov_x_y <- input$cov_icept_x_icept_y + 
            1/2*(input$cov_icept_x_slope_y + input$cov_icept_y_slope_x) +
            1/3*cov_slope_x_slope_y +
            1/12*input$mean_slope_x*input$mean_slope_y
        
        # Compute SOS
        sos <- 1 - 12*var_y*(input$mean_slope_y*var_x - input$mean_slope_x*cov_x_y)**2 /
            ((var_y*var_x - cov_x_y**2) * (12*var_x - input$mean_slope_x**2)*input$mean_slope_y**2)
        
        # And show it
        sos
    })
    
    output$lgcm_graphs <- renderPlot({
        if (input$sim_data == TRUE) {
            # Make an OpenMx model
            mx_mod <- toOpenMx(lgcm_xy)
            
            # Simulate data
            dat <- simulateData(mx_mod, input$sim_sample)
            
            dat_x <- as.data.frame(dat[, 1:(input$sim_sample/2)])
            dat_y <- as.data.frame(dat[, (input$sim_sample/2+1):input$sim_sample])
            
            names(dat_x) <- paste0("x", 1:input$timepoints)
            names(dat_y) <- paste0("y", 1:input$timepoints)
            
            dat_x_long <- gather(dat_x)
            dat_y_long <- gather(dat_y)
            
            dat_x_long$id <- rep(1:input$sim_sample, input$timepoints)
            dat_y_long$id <- rep(1:input$sim_sample, input$timepoints)
            
            names(dat_x_long)[1] <- "timepoint"
            names(dat_y_long)[1] <- "timepoint"
            
            dat_x_long$timepoint <- factor(rep(1:input$timepoints, each = input$sim_sample))
            dat_y_long$timepoint <- factor(rep(1:input$timepoints, each = input$sim_sample))
            
            dat_long <- rbind(dat_x_long, dat_y_long)
            dat_long$var <- factor(c(rep("x", nrow(dat_x_long)), 
                                     rep("y", nrow(dat_y_long))))
            
            plt_x_y <- ggplot(dat_long, aes(x = timepoint, y = value)) +
                geom_line(aes(group = id)) +
                facet_grid(var ~.)
            
            plt_x_y
        }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)