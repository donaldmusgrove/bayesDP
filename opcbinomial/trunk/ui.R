library(ggplot2)
library(knitr)
library(parallel)
library(rmarkdown)
library(shiny)
library(shinyBS)


shinyUI(pageWithSidebar(
  headerPanel("Bayesian Loss Function - Operating Characteristics"),  

  sidebarPanel(
    actionButton("button", "Submit", icon("thumbs-o-up")),
    br(),
    br(),
    wellPanel(
      h3("Simulation inputs"), 
      selectInput(inputId  = "distribution", 
                  label    = "Distribution", 
                  choices  = list("Binomial"                          = "binomial", 
                                  "Negative Binomial (not available)" = "a1", 
                                  "Gaussian (not available)"          = "a2",
                                  "Piece-wise exponential (survival)" = "a3"), 
                  selected = "binomial"),
      
      textInput(inputId = 'true_rate', 
                label   = 'True rate (use commas)', 
                value   = "0.04,0.075,0.1"),     
      
      textInput(inputId = 'weibull_scale', 
                label   = 'Weibull scale (use commas)', 
                value   = "0.01,0.05,0.5"),
      
      textInput(inputId = 'weibull_shape', 
                label   = 'Weibull shape (use commas)', 
                value   = "3"),
                
      numericInput(inputId = "sim", 
                   label   = "Number of simulations:", 
                   value   = 100, 
                   min     = 1, 
                   max     = 1000, 
                   step    = 10)        
      ),
      

    wellPanel(
      h3("Study inputs"), 
      
      numericInput(inputId = "CL", 
                   label   = "Confidence level:", 
                   value   = 0.95, 
                   min     = 0.01, 
                   max     = 0.999, 
                   step    = 0.01),
                   
      numericInput(inputId = "N", 
                   label   = "Number of current subjects:", 
                   value   = 200, 
                   min     = 1, 
                   max     = 1000, 
                   step    = 10),
                   
      numericInput(inputId = "N0", 
                   label   = "Number of historical subjects:", 
                   value   = 400, 
                   min     = 1, 
                   max     = 1000, 
                   step    = 10),
                   
      numericInput(inputId = "y0", 
                   label   = "Number of historical events:", 
                   value   = 15, 
                   min     = 1, 
                   max     = 1000, 
                   step    = 10),
                   
      numericInput(inputId = "H0_rate", 
                   label   = "Null rate:", 
                   value   = 0.1, 
                   min     = 0.01, 
                   max     = 0.999, 
                   step    = 0.01),
                   
      numericInput(inputId = "number_mcmc", 
                   label   = "Number of MCMC Simulations:", 
                   value   = 10000, 
                   min     = 50, 
                   max     = 1000000, 
                   step    = 10) 
    ),
      
    wellPanel(
      h3("Prior settings"),

      numericInput(inputId = "N0_max", 
                   label   = "Maximum effective historical sample size:", 
                   value   = 200, 
                   min     = 1, 
                   max     = 1000, 
                   step    = 10),
    
      numericInput(inputId = "alpha0", 
                   label   = "alpha - initial noninformative prior:", 
                   value   = 1, 
                   min     = 0.01, 
                   max     = 100, 
                   step    = 1),
                   
      numericInput(inputId = "beta0", 
                   label   = "beta - initial noninformative prior:", 
                   value   = 1, 
                   min     = 0.01, 
                   max     = 100, 
                   step    = 1)
    ),            
    
    wellPanel(
      h3("Download report"),
      radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
                   inline = TRUE),
      downloadButton('downloadReport')
    )
    
  ),

  mainPanel(
    tabsetPanel(
      tabPanel("Simulation Results",
               br(),
               tableOutput('sim_results')),
      tabPanel("Power Curve",
               br(),
               plotOutput("p3")),
      tabPanel("Quantile Plot",
               br(),
               plotOutput("pp3")),
      tabPanel("R Code",
               ### Currently, only displays binomial code. Want this to be dependent on user input
               includeMarkdown("www/binomialPrint.md"),
               
        ### Tooltips
        bsTooltip('distribution', "Select distribution"),
        bsTooltip('CL', "Select confidence level"),
        bsTooltip('true_rate', "Select true rate to simulate conditions under. If multiple values, separate using commas"),
        bsTooltip('H0_rate', "Select the null rate, i.e., alpha level"),
        bsTooltip('sim', "Select the number of simulations, i.e., the number of times a trial is simulated"),
        bsTooltip('N', "Select number of current subjects"),
        bsTooltip('N0', "Select number of historical subjects"),
        bsTooltip('y0', "Select number of events for historical subjects"),
        bsTooltip('NO_max', "Select maximum effective sample size the prior can receive when the loss function equals 1"),
        bsTooltip('alpha0', "Noninformative initial prior"),
        bsTooltip('beta0', "Noninformative initial prior"),
        bsTooltip('number_mcmc', "Number of simulations to estimate posterior and loss function"),
        bsTooltip('weibull_scale', "Loss function parameter controlling the location of a weibull function. If multiple values, separate using commas"),
        bsTooltip('weibull_shape', "Loss function parameter controlling the location of a weibull function. If multiple values, separate using commas"))
    )
  )

))



