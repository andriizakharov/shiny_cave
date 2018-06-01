library(shiny)
library(tidyverse)
library(OpenMx)
library(MASS)

source("input_presets.R")
source("helper_functions.R")

ui <- fluidPage(
    titlePanel("Victoria Longitudinal Study preset"),
    fluidRow(
        column(
            width = 5,
            h3("Variable X: Word Recall"),
            sliderInput(
                "mean_icept_x",
                "Mean intercept:",
                min = -1,
                max = 1,
                value = 0
            ),
            sliderInput(
                "mean_slope_x",
                "Mean slope:",
                min = -1,
                max = 1,
                value = 0
            ),
            plotOutput("lgcm_graph_x")
        ),
        column(
            width = 5,
            offset = 2,
            h3("Variable Y: Reaction Time"),
            sliderInput(
                "mean_icept_y",
                "Mean intercept:",
                min = -1,
                max = 1,
                value = 0
            ),
            sliderInput(
                "mean_slope_y",
                "Mean slope:",
                min = -1,
                max = 1,
                value = 0
            ),
            plotOutput("lgcm_graph_y")
        )
    ),
    fluidRow(actionButton("sim_data",
                          "Simulate"))
)

server <- function(input, output) {
    lgcm_x <- reactive({
        req(input$mean_icept_x, input$mean_slope_x)
        lgcm(
            timepoints = vls_RT_WRC_x_axis,
            intercept.variance = vls_RT_WRC$sigma2_Ix,
            slope.variance = vls_RT_WRC$sigma2_Sx,
            residual.variance = vls_RT_WRC$sigma2_Ex,
            intercept.slope.covariance = vls_RT_WRC$sigma_IxSx,
            intercept.mean = input$mean_icept_x,
            slope.mean = input$mean_slope_x
        )
    })
    lgcm_y <- reactive({
        req(input$mean_icept_y, input$mean_slope_y)
        lgcm(
            timepoints = vls_RT_WRC_x_axis,
            intercept.variance = vls_RT_WRC$sigma2_Iy,
            slope.variance = vls_RT_WRC$sigma2_Sy,
            residual.variance = vls_RT_WRC$sigma2_Ey,
            intercept.slope.covariance = vls_RT_WRC$sigma_IySy,
            intercept.mean = input$mean_icept_y,
            slope.mean = input$mean_slope_y
        )
    })
    lgcm_xy <- reactive({
        bivariateLgcm(
            lgcm_x(),
            lgcm_y(),
            icept.x.icept.y.covariance = vls_RT_WRC$sigma_IyIx,
            slope.x.slope.y.covariance = vls_RT_WRC$sigma_SySx
        )
    })
    
    mx_mod <- reactive({
        toOpenMx(lgcm_xy())
    })
    
    dat_long <- eventReactive(input$sim_data, {
        #req(input$sim_sample, input$timepoints)
        dat_transform1(simulateData(mx_mod(), vls_RT_WRC_n),
                       vls_RT_WRC_n,
                       vls_RT_WRC_x_axis)
    })
    
    output$lgcm_graph_x <- renderPlot({
        ggplot(dat_long()[[1]], aes(x = timepoint, y = value)) +
            geom_line(aes(group = id))
    })
    
    output$lgcm_graph_y <- renderPlot({
        ggplot(dat_long()[[2]], aes(x = timepoint, y = value)) +
            geom_line(aes(group = id))
    })
    
}





### Run the app
shinyApp(ui = ui, server = server)