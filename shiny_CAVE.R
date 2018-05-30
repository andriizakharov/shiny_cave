library(semper)
source("bivariateLgcm.R")
source("toOpenMx.R")
#source("input_presets.R")
source("helper_functions.R")

library(shiny)
library(ggplot2)
library(tidyverse)


ui <- fluidPage(
    
    # Application title
    titlePanel("Bivariate LGCM simulation"),
    
    # Sidebar
    sidebarLayout(
        sidebarPanel(
            h2("Presets:"),
            
            selectInput("preset",
                        "Choose a preset:",
                        choices = list("",
                            "custom"
                        )),
            
            conditionalPanel("input.preset == 'custom'",
                h2("Variable X:"),
            
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
                
                h2("Variable Y:"),
                
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
                
                h2("Cross-variable:"),
                
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
                
                h2("Other:"),
                
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
                             "Simulate!"),
                
                # Set seed for sampling
                numericInput("seed",
                             "Set random seed for sampling",
                             min = -100,
                             max = 100,
                             value = 42),
                
                # Button to draw a cross-sectional sample
                actionButton("sample",
                             "Sample cross-section")
            )
        ),
        
        # Show 
        mainPanel(
            titlePanel("SOS:"),
            h2(htmlOutput("sos_val")),
            plotOutput("lgcm_graphs"),
            plotOutput("scatter_plots")
        )
    )
)


server <- function(input, output) {
    
    ### Reactive calculations
    
    # Define the models given all input parameters
    lgcm_x <- reactive({ 
        req(input$timepoints, input$var_icept_x, input$var_slope_x,
            input$error_var_x, input$cov_icept_x_slope_x, input$mean_icept_x,
            input$mean_slope_x)
        lgcm(
            timepoints = 1:input$timepoints,
            intercept.variance = input$var_icept_x,
            slope.variance = input$var_slope_x,
            residual.variance = input$error_var_x,
            intercept.slope.covariance = input$cov_icept_x_slope_x,
            intercept.mean = input$mean_icept_x,
            slope.mean = input$mean_slope_x
    ) })
    lgcm_y <- reactive({ 
        req(input$timepoints, input$var_icept_y, input$var_slope_y,
            input$error_var_y, input$cov_icept_y_slope_y, input$mean_icept_y,
            input$mean_slope_y)
        lgcm(
            timepoints = 1:input$timepoints,
            intercept.variance = input$var_icept_y,
            slope.variance = input$var_slope_y,
            residual.variance = input$error_var_y,
            intercept.slope.covariance = input$cov_icept_y_slope_y,
            intercept.mean = input$mean_icept_y,
            slope.mean = input$mean_slope_y
    ) })
    
    cov_slope_x_slope_y <- reactive({ 
        req(input$cor_slope_x_slope_y, input$var_slope_x, input$var_slope_y)
        input$cor_slope_x_slope_y * sqrt(input$var_slope_x * input$var_slope_y) 
    })
    
    lgcm_xy <- reactive({ 
        req(input$cov_icept_x_icept_y)
        bivariateLgcm(
            lgcm_x(), lgcm_y(), 
            icept.x.icept.y.covariance = input$cov_icept_x_icept_y, 
            slope.x.slope.y.covariance = cov_slope_x_slope_y()
        ) })
    
    var_x <- reactive({ 
        req(input$var_icept_x, input$var_slope_x, input$cov_icept_x_slope_x,
            input$mean_slope_x, input$error_var_x)
        input$var_icept_x + 1/3*input$var_slope_x + input$cov_icept_x_slope_x + 
        1/12*input$mean_slope_x**2 + input$error_var_x })
    var_y <- reactive({ 
        req(input$var_icept_y, input$var_slope_y, input$cov_icept_y_slope_y,
            input$mean_slope_y, input$error_var_y)
        input$var_icept_y + 1/3*input$var_slope_y + input$cov_icept_y_slope_y + 
        1/12*input$mean_slope_y**2 + input$error_var_y })
    
    cov_x_y <- reactive({ 
        req(input$cov_icept_x_icept_y, input$cov_icept_x_slope_y, 
            input$cov_icept_y_slope_x, input$mean_slope_x, input$mean_slope_y)
        input$cov_icept_x_icept_y + 
        1/2*(input$cov_icept_x_slope_y + input$cov_icept_y_slope_x) +
        1/3*cov_slope_x_slope_y() +
        1/12*input$mean_slope_x*input$mean_slope_y })
    
    # Compute SOS
    sos <- reactive({ 
        req(input$mean_slope_y, input$mean_slope_x)
        1 - 12*var_y()*(input$mean_slope_y*var_x() - input$mean_slope_x*cov_x_y())**2 /
        ((var_y()*var_x() - cov_x_y()**2) * (12*var_x() - input$mean_slope_x**2)*
             input$mean_slope_y**2) })
    
    # Make an OpenMx model
    mx_mod <- reactive({ toOpenMx(lgcm_xy()) })
    
    # Simulate data if button is clicked
    dat_long <- eventReactive(input$sim_data, {
        req(input$sim_sample, input$timepoints)
        dat_transform(simulateData(mx_mod(), input$sim_sample), 
                      input$sim_sample, input$timepoints)
    })
    
    dat_one_per_person <- eventReactive(input$sample, {
        dat_sample(dat_long(), input$sim_sample, seed = input$seed)
    })
    
    output$sos_val <- renderText( sos() )
    
    output$lgcm_graphs <- renderPlot({
        
        ggplot(dat_long(), aes(x = timepoint, y = value)) +
            geom_line(aes(group = id)) +
            facet_grid(var ~.)
        
    })
    
    output$scatter_plots <- renderPlot({
        
        ggplot(dat_one_per_person(), aes(x = timepoint, y = value)) +
            geom_point() +
            facet_grid(var ~.)
    })
}




# Run the application 
shinyApp(ui = ui, server = server)