library(shiny)


shinyUI(pageWithSidebar(
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
               includeMarkdown("gaussianPrint.md")
      )               
    )
  )
))
