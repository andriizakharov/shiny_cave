library(shiny)
library(ggplot2)
#library(OpenMx)
#library(MASS)

source("input_presets.R")
source("helper_functions.R")

ui <- fluidPage(
    titlePanel("Slope-slope correlation vs. Shared-over-simple effects"),
    fluidRow(
        column(width = 6,
            h2("Presets:"),
            selectInput("preset",
                        "Choose a preset:",
                        choices = list("LOGH_2011_A", "LOGH_2011_B", "LOGH_2011_C",
                                       "custom", "vls_RT_WRC", "acad_A_CD", "elsa_AF_PM"
                        ))
        )
    ),
    
    conditionalPanel("input.preset == 'vls_RT_WRC' || 
                     input.preset == 'acad_A_CD' ||
                     input.preset == 'elsa_AF_PM'",
        numericInput("mean_slope_x",
                     "X: slope mean",
                     value = 1),
        numericInput("mean_slope_y",
                     "Y: slope mean",
                     value = 1)
    ),
    
    conditionalPanel("input.preset == 'custom'",
        fluidRow(
            column(width = 6,
                h2("Variable X:"),
        
                # Input params for variable X
                # numericInput("mean_icept_x",
                #           "X: intercept mean",
                #           value = 0),
                numericInput("mean_slope_x",
                          "X: slope mean",
                          value = 1),
                numericInput("sigma2_Ix",
                          "X: intercept variance",
                          value = 0),
                numericInput("sigma2_Sx",
                          "X: slope variance",
                          value = 0),
                numericInput("sigma2_Ex",
                          "X: error variance",
                          value = 0),
                numericInput("sigma_IxSx",
                          "X: intercept-slope covariance",
                          value = 0)
            ),
        
            column(width = 6,
                h2("Variable Y:"),
        
                # Input params for variable Y
                # numericInput("mean_icept_y",
                #           "Y: intercept mean",
                #           value = 0),
                numericInput("mean_slope_y",
                          "Y: slope mean",
                          value = 1),
                numericInput("sigma2_Iy",
                          "Y: intercept variance",
                          value = 0),
                numericInput("sigma2_Sy",
                          "Y: slope variance",
                          value = 0),
                numericInput("sigma2_Ey",
                          "Y: error variance",
                          value = 0),
                numericInput("sigma_IySy",
                          "Y: intercept-slope covariance",
                          value = 0)
            )
        ),
        
        fluidRow(
            column(width = 6,
                h2("Cross-variable:"),
        
                # Input intervariable params
                numericInput("sigma_IyIx",
                          "XY: intercept covariance",
                          value = 0),
                numericInput("sigma_SyIx",
                          "XY: covariance intercept X - slope Y",
                          value = 0),
                numericInput("sigma_IySx",
                          "XY: covariance intercept Y - slope X",
                          value = 0)
                
                # sliderInput("cor_slope_x_slope_y",
                #          "XY: slope-slope correlation",
                #          min = 0,
                #          max = 1,
                #          value = 0.5)
            )

            # column(width = 5,
            #     h2("Other:"),
            # 
            #     # How many timepoints should the models simulate?
            #     numericInput("timepoints",
            #               "Number of timepoints per model",
            #               min = 0,
            #               max = 100,
            #               value = 10),
            #     
            #     # How many simulated trajectories?
            #     numericInput("sim_sample",
            #               "Number of simulated trajectories",
            #               min = 0,
            #               max = 100,
            #               value = 20)
            # )
        )
        
        # Button to simulate data
        #actionButton("sim_data",
        #          "Simulate!")
    ),
    fluidRow(
        column(width = 12,
               tableOutput("params")
            
        )
    ),
    fluidRow(plotOutput("sos_plot")),
    fluidRow(column(width = 6,
                    actionButton("sos",
                          "SOS!"))
    )
    
)

#     fluidRow(
#         column(
#             width = 5,
#             h3("Variable X: Word Recall"),
#             sliderInput(
#                 "mean_icept_x",
#                 "Mean intercept:",
#                 min = 0,
#                 max = 100,
#                 value = 0
#             ),
#             sliderInput(
#                 "mean_slope_x",
#                 "Mean slope:",
#                 min = -20,
#                 max = 20,
#                 value = 75
#             ),
#             plotOutput("lgcm_graph_x")
#         ),
#         column(
#             width = 5,
#             offset = 2,
#             h3("Variable Y: Reaction Time"),
#             sliderInput(
#                 "mean_icept_y",
#                 "Mean intercept:",
#                 min = 0,
#                 max = 100,
#                 value = 0
#             ),
#             sliderInput(
#                 "mean_slope_y",
#                 "Mean slope:",
#                 min = -2000,
#                 max = 2000,
#                 value = 70
#             ),
#             plotOutput("lgcm_graph_y")
#         )
#     ),
#     fluidRow(actionButton("sim_data",
#                           "Simulate")),
#     fluidRow(plotOutput("sos_plot")),
#     fluidRow(actionButton("sos",
#                           "SOS!"))
# )

server <- function(input, output) {
    # lgcm_x <- reactive({
    #     req(input$mean_icept_x, input$mean_slope_x)
    #     lgcm(
    #         timepoints = vls_RT_WRC_x_axis,
    #         intercept.variance = vls_RT_WRC$sigma2_Ix,
    #         slope.variance = vls_RT_WRC$sigma2_Sx,
    #         residual.variance = vls_RT_WRC$sigma2_Ex,
    #         intercept.slope.covariance = vls_RT_WRC$sigma_IxSx,
    #         intercept.mean = input$mean_icept_x,
    #         slope.mean = input$mean_slope_x
    #     )
    # })
    # lgcm_y <- reactive({
    #     req(input$mean_icept_y, input$mean_slope_y)
    #     lgcm(
    #         timepoints = vls_RT_WRC_x_axis,
    #         intercept.variance = vls_RT_WRC$sigma2_Iy,
    #         slope.variance = vls_RT_WRC$sigma2_Sy,
    #         residual.variance = vls_RT_WRC$sigma2_Ey,
    #         intercept.slope.covariance = vls_RT_WRC$sigma_IySy,
    #         intercept.mean = input$mean_icept_y,
    #         slope.mean = input$mean_slope_y
    #     )
    # })
    # lgcm_xy <- reactive({
    #     bivariateLgcm(
    #         lgcm_x(),
    #         lgcm_y(),
    #         icept.x.icept.y.covariance = vls_RT_WRC$sigma_IyIx,
    #         slope.x.slope.y.covariance = vls_RT_WRC$sigma_SySx
    #     )
    # })
    # 
    # mx_mod <- reactive({
    #     toOpenMx(lgcm_xy())
    # })
    # 
    # dat_long <- eventReactive(input$sim_data, {
    #     #req(input$sim_sample, input$timepoints)
    #     dat_transform1(simulateData(mx_mod(), vls_RT_WRC_n),
    #                    vls_RT_WRC_n,
    #                    vls_RT_WRC_x_axis)
    # })
    # 
    # output$lgcm_graph_x <- renderPlot({
    #     ggplot(dat_long()[[1]], aes(x = timepoint, y = value)) +
    #         geom_line(aes(group = id)) +
    #         xlab("Age")
    # })
    # 
    # output$lgcm_graph_y <- renderPlot({
    #     ggplot(dat_long()[[2]], aes(x = timepoint, y = value)) +
    #         geom_line(aes(group = id)) +
    #         xlab("Age")
    # })
    
    # sos_cor_df <- eventReactive(input$sos, {
    #     req(input$mean_slope_x, input$mean_slope_y)
    #     sos_by_cor(input$mean_slope_x, vls_RT_WRC$sigma2_Ix, vls_RT_WRC$sigma2_Sx,
    #                vls_RT_WRC$sigma2_Ex,
    #                input$mean_slope_y, vls_RT_WRC$sigma2_Iy, vls_RT_WRC$sigma2_Sy,
    #                vls_RT_WRC$sigma2_Ey,
    #                vls_RT_WRC$sigma_IxSx, vls_RT_WRC$sigma_IySy, vls_RT_WRC$sigma_SyIx,
    #                vls_RT_WRC$sigma_IySx, vls_RT_WRC$sigma_IyIx)
    # })
    sos_cor_df <- eventReactive(input$sos, {
        if(input$preset == "LOGH_2011_A") {
            sos_by_cor(LOGH_2011_A$mean_slope_x, LOGH_2011_A$sigma2_Ix, LOGH_2011_A$sigma2_Sx,
                       LOGH_2011_A$sigma2_Ex,
                       LOGH_2011_A$mean_slope_y, LOGH_2011_A$sigma2_Iy, LOGH_2011_A$sigma2_Sy,
                       LOGH_2011_A$sigma2_Ey,
                       LOGH_2011_A$sigma_IxSx, LOGH_2011_A$sigma_IySy, LOGH_2011_A$sigma_SyIx,
                       LOGH_2011_A$sigma_IySx, LOGH_2011_A$sigma_IyIx)

        } else if(input$preset == "LOGH_2011_B"){
            sos_by_cor(LOGH_2011_B$mean_slope_x, LOGH_2011_B$sigma2_Ix, LOGH_2011_B$sigma2_Sx,
                       LOGH_2011_B$sigma2_Ex,
                       LOGH_2011_B$mean_slope_y, LOGH_2011_B$sigma2_Iy, LOGH_2011_B$sigma2_Sy,
                       LOGH_2011_B$sigma2_Ey,
                       LOGH_2011_B$sigma_IxSx, LOGH_2011_B$sigma_IySy, LOGH_2011_B$sigma_SyIx,
                       LOGH_2011_B$sigma_IySx, LOGH_2011_B$sigma_IyIx)
        
        } else if(input$preset == "LOGH_2011_C"){
            sos_by_cor(LOGH_2011_C$mean_slope_x, LOGH_2011_C$sigma2_Ix, LOGH_2011_C$sigma2_Sx,
                       LOGH_2011_C$sigma2_Ex,
                       LOGH_2011_C$mean_slope_y, LOGH_2011_C$sigma2_Iy, LOGH_2011_C$sigma2_Sy,
                       LOGH_2011_C$sigma2_Ey,
                       LOGH_2011_C$sigma_IxSx, LOGH_2011_C$sigma_IySy, LOGH_2011_C$sigma_SyIx,
                       LOGH_2011_C$sigma_IySx, LOGH_2011_C$sigma_IyIx)
        
        } else if(input$preset == "custom"){
            sos_by_cor(input$mean_slope_x, input$sigma2_Ix, input$sigma2_Sx,
                       input$sigma2_Ex,
                       input$mean_slope_y, input$sigma2_Iy, input$sigma2_Sy,
                       input$sigma2_Ey,
                       input$sigma_IxSx, input$sigma_IySy, input$sigma_SyIx,
                       input$sigma_IySx, input$sigma_IyIx)
        
        } else if(input$preset == "vls_RT_WRC"){
            sos_by_cor(input$mean_slope_x, vls_RT_WRC$sigma2_Ix, vls_RT_WRC$sigma2_Sx,
                       vls_RT_WRC$sigma2_Ex,
                       input$mean_slope_y, vls_RT_WRC$sigma2_Iy, vls_RT_WRC$sigma2_Sy,
                       vls_RT_WRC$sigma2_Ey,
                       vls_RT_WRC$sigma_IxSx, vls_RT_WRC$sigma_IySy, vls_RT_WRC$sigma_SyIx,
                       vls_RT_WRC$sigma_IySx, vls_RT_WRC$sigma_IyIx)
        
        } else if(input$preset == "acad_A_CD"){
            sos_by_cor(input$mean_slope_x, acad_A_CD$sigma2_Ix, acad_A_CD$sigma2_Sx,
                       acad_A_CD$sigma2_Ex,
                       input$mean_slope_y, acad_A_CD$sigma2_Iy, acad_A_CD$sigma2_Sy,
                       acad_A_CD$sigma2_Ey,
                       acad_A_CD$sigma_IxSx, acad_A_CD$sigma_IySy, acad_A_CD$sigma_SyIx,
                       acad_A_CD$sigma_IySx, acad_A_CD$sigma_IyIx)
        
        } else if(input$preset == "elsa_AF_PM"){
            sos_by_cor(input$mean_slope_x, elsa_AF_PM$sigma2_Ix, elsa_AF_PM$sigma2_Sx,
                       elsa_AF_PM$sigma2_Ex,
                       input$mean_slope_y, elsa_AF_PM$sigma2_Iy, elsa_AF_PM$sigma2_Sy,
                       elsa_AF_PM$sigma2_Ey,
                       elsa_AF_PM$sigma_IxSx, elsa_AF_PM$sigma_IySy, elsa_AF_PM$sigma_SyIx,
                       elsa_AF_PM$sigma_IySx, elsa_AF_PM$sigma_IyIx)
        }
        
    })
    
    output$sos_plot <- renderPlot({
        ggplot(sos_cor_df(), aes(x = cor, y = sos)) +
            geom_line() +
            scale_y_continuous("Shared-over-simple effects",
                               limits = c(0,1)) +
            labs(title = "Shared-over-simple (SOS) effects by slope-slope correlation 
                  of the variables X and Y", x = "Slope-slope correlation X-Y") 
    })
    
    output$params <- renderTable({
        get(input$preset)
    })
    
}





### Run the app
shinyApp(ui = ui, server = server)