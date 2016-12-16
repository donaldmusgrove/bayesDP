######################################################################
# Create the first quadrant class
#
# This is used to represent a coordinate in the first quadrant.
rtcnormal <- setClass(
  # Set the name for the class
  "rtcnormal",

  # Define the slots
  slots = c(
    ui = "character",
    server = "function"
  )#,

  # Set the default values for the slots. (optional)
  #prototype=list(
  #  x = 0.0,
  #  y = 0.0
  #),

  # Make a function that can test to see if the data is consistent.
  # This is not called if you have an initialize function defined!
  #validity=function(object)
  #{
  #  if((object@x < 0) || (object@y < 0)) {
  #    return("A negative number for one of the coordinates was given.")
  #  }
  #  return(TRUE)
  #}
)

setGeneric(name="shine",
           def=function(x)
           {
             standardGeneric("shine")
           }
)

setMethod(f="shine",
          signature="rtcnormal ",
          definition=function(x)
          {

################################################################################
# UI
################################################################################

ui <- shinyUI(pageWithSidebar(
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
            includeMarkdown("./inst/www/gaussianPrint.md"))
  ))
))


################################################################################
# Server Code
################################################################################

server <- shinyServer(function(input, output, session){
  source("./R/post_estimation_normal_functions.R")
  source("./R/Normal_RCT_shiny.R")


  est <- reactive({
    mu_t <- mu_posterior(
      mu            = input$mu,
      sigma2        = input$sigma2,
      N             = input$N,             #Number of  current subjects
      mu0           = input$mu0,
      sigma02       = input$sigma02 ,      #Number of events observed  historical  data sets
      N0            = input$N0,            #Number of historical subjects
      N0_max        = input$N0_max,        #Maximum effective sample size the prior can receive when the loss function equals 1
      number_mcmc   = input$number_mcmc,   #Number of simulations to estimate posterior and loss function
      weibull_scale = input$weibull_scale, #Loss function parameter controlling the location of a weibull function
      weibull_shape = input$weibull_shape, #Loss function parameter controlling the location of a weibull function
      two_side      = input$two_side)

    mu_c <- mu_posterior(
      mu            = input$mu_c,
      sigma2        = input$sigma2_c,
      N             = input$N_c,             #Number of  current subjects
      mu0           = input$mu0_c,
      sigma02       = input$sigma02_c,       #Number of events observed  historical  data sets
      N0            = input$N0_c,            #Number of historical subjects
      N0_max        = input$N0_max_c,        #Max effective sample size the prior can receive when the loss function equals 1
      number_mcmc   = input$number_mcmc,     #Number of simulations to estimate posterior and loss function
      weibull_scale = input$weibull_scale_c, #Loss function parameter controlling the location of a weibull function
      weibull_shape = input$weibull_shape_c, #Loss function parameter controlling the location of a weibull function
      two_side      = input$two_side)

    return(list(mu_t = mu_t,
                mu_c = mu_c))
  })


  f1 <- reactive({
    final(posterior_test    = est()$mu_t,
          posterior_control = est()$mu_c)
  })


  res1 <- reactive({
    results(f                 = f1(),
            posterior_test    = est()$mu_t,
            posterior_control = est()$mu_c,
            two_side          = input$two_side,
            inequality        = input$inequality,
            N0_t              = input$N0,
            N0_c              = input$N0_c,
            delta             = input$delta)
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

  output$hypothesis1<-renderPrint({
    cat(res1()$hypothesis)
  })

  output$prior_for_test_group1<-renderPrint({
    res1()$prior_for_test_group
  })

  output$prior_for_control_group1<-renderPrint({
    res1()$prior_for_control_group
  })
})

pane_runApp(ui,server)
          }
)

