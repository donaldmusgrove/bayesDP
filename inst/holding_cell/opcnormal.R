library(shiny)
library(shinyFiles)
library(shinythemes)
library(shinyBS)

library(survival)
library(ggplot2)
library(MCMCpack)
library(knitr)
library(rmarkdown)
library(parallel)

.runApp <- function(app, ...) {
  ## selectively use the RStudio viewer pane (if available)
  viewer <- getOption("viewer")
  if (!is.null(viewer) && is.function(viewer)) {
    runApp(app, launch.browser = viewer, ...)
  } else {
    runApp(app, ...)
  }
}

################################################################################
# UI
################################################################################

options(shiny.launch.browser = .rs.invokeShinyWindowViewer)


ui <- shinyUI(pageWithSidebar(
  headerPanel("NormalOPC"),
  sidebarPanel(
    wellPanel(
      h4("loss fuction info"),
      numericInput("N0_max", "N0_max", value=10,min =1, max =100,step = 1),
      checkboxInput('two_side', 'two_side', FALSE),
      numericInput("weibull_scale", "weibull_scale", value=.1,min =1, max =100,step = 1),
      numericInput("weibull_shape", "weibull_shape", value=2,min =1, max =100,step = 1)
    ),

    wellPanel(
      h4("Current data info"),
      numericInput("mu", "mu", value=10,min =1, max =100,step = 1),
      numericInput("sigma2", "sigma2", value=1,min =1, max =100,step = 1),
      numericInput("N", "N", value=10,min =1, max =100,step = 1)
    ),

    wellPanel(
      h4("hist data info"),
      numericInput("mu0", "mu0", value=10,min =1, max =100,step = 1),
      numericInput("sigma02", "sigma02", value=1,min =1, max =100,step = 1),
      numericInput("N0", "N0", value=10,min =1, max =100,step = 1)
    ),

    numericInput("number_mcmc", "number_mcmc", value=1000,min =1, max =100,step = 1),

    tags$hr(),
    tags$style(type='text/css', "#downloadData { vertical-align: top;}")
  ),

  mainPanel(
    tabsetPanel(

      tabPanel("hypothesis_test",
               wellPanel(
                 h4("Hypothesis Test info"),
                 numericInput("H0", "H0 Value", value=10,min =1, max =100,step = 1),
                 selectInput("inequality",
                             label = h3("Inequality of alternate hypothesis (H_a)"),
                             choices = list(">" =">", "<" = "<"),
                             selected = 1),
                 verbatimTextOutput('hypothesis1')
               ),
               h4("Strength of prior data:",verbatimTextOutput('prior_for_test_group1')),plotOutput("lossfun_plot2")
      ),

      tabPanel("plot",
               plotOutput("post_typeplot1"),
               plotOutput("lossfun_plot1"),
               plotOutput("densityplot1")
      ),

      tabPanel("R Code",
               includeMarkdown("./inst/www/gaussianPrint.md")
      )
    )
  )
))

################################################################################
# Server Code
################################################################################

server <- shinyServer(function(input, output, session){

  source("./R/post_estimation_normal_functions.R")
  source("./R/Normal_OPC_shiny.R")

  est <- reactive({
    estPost <- mu_posterior(mu            = input$mu,
                            sigma2        = input$sigma2,
                            N             = input$N,             #Number of  current subjects
                            mu0           = input$mu0,
                            sigma02       = input$sigma02,
                            N0            = input$N0,            #Number of historical subjects
                            N0_max        = input$N0_max,        #Maximum effective prior sample size
                            number_mcmc   = input$number_mcmc,   #Number of posterior simulations
                            weibull_scale = input$weibull_scale, #Loss function: Weibull location scale
                            weibull_shape = input$weibull_shape, #Loss function: Weibull location shape
                            two_side      = input$two_side)
    estPost
  })

  f1 <- reactive({
    final1(posterior_test = est())
  })

  res1 <- reactive({
    results(f              = f1(),
            posterior_test = est(),
            H0             = input$H0,
            two_side       = input$two_side,
            inequality     = input$inequality)
  })

  output$post_typeplot1 <- renderPlot({
    print(res1()$post_typeplot)
  })

  output$densityplot1 <- renderPlot({
    print(res1()$densityplot)
  })

  output$lossfun_plot1 <- renderPlot({
    print(res1()$lossfun_plot)
  })

  output$lossfun_plot2 <- renderPlot({
    print(res1()$lossfun_plot)
  })

  output$hypothesis1 <- renderPrint({
    cat(res1()$hypothesis)
  })

  output$prior_for_test_group1 <- renderPrint({
    res1()$prior_for_test_group
  })
})

################################################################################
# Finally Run Everything
################################################################################

.runApp(shinyApp(ui, server))
#runApp(shinyApp(ui, server),
#       port=as.numeric(commandArgs(TRUE)[2]),
#       launch.browser=FALSE)

