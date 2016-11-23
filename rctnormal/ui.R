library(shiny)


shinyUI(pageWithSidebar(
   # Application title
  headerPanel("NormalRCT"),  
  sidebarPanel( 
    
    wellPanel(h2('Loss Function Inputs:') ,
      wellPanel( 
        h4("Test group info"),
        numericInput("N0_max", "N0_max", value=10,min =1, max =100,step = 1),
        numericInput("weibull_scale", "weibull_scale", value=.1,min =1, max =100,step = 1),
        numericInput("weibull_shape", "weibull_shape", value=2,min =1, max =100,step = 1)
      ),
      
      wellPanel( 
        h4("Control group info"),
        numericInput("N0_max_c", "N0_max", value=10,min =1, max =100,step = 1),
        numericInput("weibull_scale_c", "weibull_scale", value=.1,min =1, max =100,step = 1),
        numericInput("weibull_shape_c", "weibull_shape", value=2,min =1, max =100,step = 1)
      ), 
      
      checkboxInput('two_side', 'two side loss function', FALSE)
    ),
    
    wellPanel( 
      h2("test group info"),
      wellPanel( 
        h4("Current data info"),
           # checkboxInput('two_side', 'two_side', FALSE),
        numericInput("mu", "mu", value=10,min =1, max =100,step = 1),
        numericInput("sigma2", "sigma2", value=1,min =1, max =100,step = 1),
        numericInput("N", "N", value=10,min =1, max =100,step = 1)
      ),
      wellPanel( 
        h4("Hist data info"),
        numericInput("mu0", "mu0", value=10,min =1, max =100,step = 1),
        numericInput("sigma02", "sigma02", value=1,min =1, max =100,step = 1),
        numericInput("N0", "N0", value=10,min =1, max =100,step = 1)
      )
    ),
    
    wellPanel( 
      h2("Control group info"),
      wellPanel( 
        h4("Current data info"),
           # checkboxInput('two_side', 'two_side', FALSE),
        numericInput("mu_c", "mu", value=10,min =1, max =100,step = 1),
        numericInput("sigma2_c", "sigma2", value=1,min =1, max =100,step = 1),
        numericInput("N_c", "N", value=10,min =1, max =100,step = 1)
      ),
      wellPanel( 
        h4("Hist data info"),
        numericInput("mu0_c", "mu0", value=10,min =1, max =100,step = 1),
        numericInput("sigma02_c", "sigma02", value=1,min =1, max =100,step = 1),
        numericInput("N0_c", "N0", value=10,min =1, max =100,step = 1)
      )
    ),

    numericInput("number_mcmc", "number_mcmc", value=1000,min =1, max =100,step = 1)
  ),


  mainPanel(
    tabsetPanel(
      tabPanel("hypothesis_test", 
        wellPanel( 
          h4("Hypothesis Test info"),
          numericInput("delta", "Non-inferiority zone value (delta) ", value=0,min =-100, max =100,step = 1), 
          selectInput("inequality", label = h3("Inequality of alternate hypothesis (H_a)"), 
                      choices = list(">" =">", "<" = "<"), 
                      selected = 1),
          verbatimTextOutput('hypothesis1')
        ),

        h2("Strength of prior data:"), 
        h5("test group"),
        verbatimTextOutput('prior_for_test_group1'), 
        h5("control group"),
        verbatimTextOutput('prior_for_control_group1'),
        plotOutput("lossfun_plot2")
      ),   
      tabPanel("plot",
               plotOutput("post_typeplot1"),
               plotOutput("lossfun_plot1"),
               plotOutput("densityplot1")
      ),
      tabPanel("R Code",
               includeMarkdown("gaussianPrint.md"))
))
   ))
