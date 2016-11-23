library(ggplot2)
library(knitr)
library(parallel)
library(rmarkdown)
library(shiny)
library(shinyBS)


shinyUI(pageWithSidebar(
  headerPanel("Bayesian Loss Function - RCT"),  

  sidebarPanel(
    br(),
    wellPanel(
      h3("Test inputs"), 
      numericInput(inputId = "y_test", 
                   label   = "Number of observed events (current data):", 
                   value   = 10, 
                   min     = 1, 
                   max     = 1000, 
                   step    = 1),   

      numericInput(inputId = "N_test", 
                   label   = "Number of current subjects:", 
                   value   = 400, 
                   min     = 1, 
                   max     = 1000, 
                   step    = 1), 
                   
       numericInput(inputId = "y0_test", 
                   label   = "Number of observed events (historical data):", 
                   value   = 1, 
                   min     = 1, 
                   max     = 100, 
                   step    = 1),                   

      numericInput(inputId = "N0_test", 
                   label   = "Number of historical subjects:", 
                   value   = 1000, 
                   min     = 1, 
                   max     = 10000, 
                   step    = 1),   

      numericInput(inputId = "N0_max_test", 
                   label   = "Max effective sample size:", 
                   value   = 200, 
                   min     = 1, 
                   max     = 1000, 
                   step    = 1), 
                   
       numericInput(inputId = "alpha0_test", 
                   label   = "alpha - noninformative initial prior:", 
                   value   = 1, 
                   min     = 0.01, 
                   max     = 100, 
                   step    = 0.1), 
 
       numericInput(inputId = "beta0_test", 
                   label   = "beta - noninformative initial prior:", 
                   value   = 1, 
                   min     = 0.01, 
                   max     = 100, 
                   step    = 0.1),                   

      numericInput(inputId = "weibull_scale_test", 
                   label   = "Weibull scale:", 
                   value   = 0.2, 
                   min     = 0.01, 
                   max     = 100, 
                   step    = 0.1),   

      numericInput(inputId = "weibull_shape_test", 
                   label   = "Weibull shape:", 
                   value   = 2, 
                   min     = 1, 
                   max     = 1000, 
                   step    = 1)
    ),
      
    wellPanel(
      h3("Control inputs"), 
      numericInput(inputId = "y_control", 
                   label   = "Number of observed events (current data):", 
                   value   = 1, 
                   min     = 1, 
                   max     = 1000, 
                   step    = 1),   

      numericInput(inputId = "N_control", 
                   label   = "Number of current subjects:", 
                   value   = 400, 
                   min     = 1, 
                   max     = 1000, 
                   step    = 1), 
                   
       numericInput(inputId = "y0_control", 
                   label   = "Number of observed events (historical data):", 
                   value   = 10, 
                   min     = 1, 
                   max     = 100, 
                   step    = 1),                   

      numericInput(inputId = "N0_control", 
                   label   = "Number of historical subjects:", 
                   value   = 10, 
                   min     = 1, 
                   max     = 10000, 
                   step    = 1),   

      numericInput(inputId = "N0_max_control", 
                   label   = "Max effective sample size:", 
                   value   = 200, 
                   min     = 1, 
                   max     = 1000, 
                   step    = 1), 
                   
       numericInput(inputId = "alpha0_control", 
                   label   = "alpha - noninformative initial prior:", 
                   value   = 1, 
                   min     = 0.01, 
                   max     = 100, 
                   step    = 0.1), 
 
       numericInput(inputId = "beta0_control", 
                   label   = "beta - noninformative initial prior:", 
                   value   = 1, 
                   min     = 0.01, 
                   max     = 100, 
                   step    = 0.1),                   

      numericInput(inputId = "weibull_scale_control", 
                   label   = "Weibull scale:", 
                   value   = 0.05, 
                   min     = 0.01, 
                   max     = 100, 
                   step    = 0.1),   

      numericInput(inputId = "weibull_shape_control", 
                   label   = "Weibull shape:", 
                   value   = 2, 
                   min     = 1, 
                   max     = 1000, 
                   step    = 1)
    ),
    
    checkboxInput('two_side', 'two side loss function', FALSE),
    
    numericInput(inputId = "number_mcmc", 
                 label   = "Number of MCMC simulations:", 
                 value   = 10000, 
                 min     = 100, 
                 max     = 50000, 
                 step    = 1)#,

    #wellPanel(
    #  h3("Download report"),
    #  radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
    #               inline = TRUE),
    #  downloadButton('downloadReport')
    #)
    
  ),

  mainPanel(
    tabsetPanel(
      tabPanel("Results",
               br(),
               uiOutput("text_test"),
               br(),
               uiOutput("text_control"),
               br(),
               uiOutput("hypothesis_text")),
      tabPanel("Densities",
               br(),
               h3("Prior vs Posterior"),
               plotOutput("plotP"),
               br(),
               h3("Test vs Control"),
               plotOutput("plotP1")),
      tabPanel("Histogram",
               br(),
               plotOutput("plotHist")),
      tabPanel("Loss Function",
               br(),
               plotOutput("plotP3")),
      tabPanel("R Code",
               ### Currently, only displays binomial code. Want this to be dependent on user input
               includeMarkdown("www/binomialPrint.md"),
               
        ### Tooltips
        bsTooltip('y_test', "Number of observed events in the current test group"),
        bsTooltip('N_test', "Number of current subjects in the test group"),
        bsTooltip('y0_test', "Number of observed events in the historical test group"),
        bsTooltip('N0_test', "Number of historical subjects in the test group"),
        bsTooltip('N0_max_test', "Maximum effective sample size"),
        bsTooltip('alpha0_test', "Noninformative initial prior"),
        bsTooltip('beta0_test', "Noninformative initial prior"),
        bsTooltip('weibull_scale_test', "Weibull scale parameter"),
        bsTooltip('weibull_shape_test', "Weibull shape parameter"),
        bsTooltip('y_control', "Number of observed events in the current control group"),
        bsTooltip('N_control', "Number of current subjects in the control group"),
        bsTooltip('y0_control', "Number of observed events in the historical control group"),
        bsTooltip('N0_test', "Number of historical subjects in the control group"),
        bsTooltip('N0_max_control', "Maximum effective sample size"),
        bsTooltip('alpha0_control', "Noninformative initial prior"),
        bsTooltip('beta0_control', "Noninformative initial prior"),
        bsTooltip('weibull_scale_control', "Weibull scale parameter"),
        bsTooltip('weibull_shape_control', "Weibull shape parameter"),
        bsTooltip('number_mcmc', "Number of simulations to estimate posterior and loss functions"))
    )
  )

))



