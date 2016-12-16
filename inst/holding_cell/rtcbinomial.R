################################################################################
# UI
################################################################################

ui <- shinyUI(pageWithSidebar(
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
               includeMarkdown("./rctbinomial/trunk/www/binomialPrint.md"),

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

################################################################################
# Server Code
################################################################################

server <- shinyServer(function(input, output, session){


  ### Load external R functions based on distribution choice
  #observeEvent(input$distribution,{
  #  if(input$distribution %in% c("binomial","a1","a2","a3")){
  source("./rctbinomial/trunk/www/post_estimation_binomial_functions.R")
  #  }
  #})



  ### Fit test data
  posterior_test <- reactive({
    Binomial_posterior(y             = input$y_test,
                       N             = input$N_test,
                       y0            = input$y0_test,
                       N0            = input$N0_test,
                       N0_max        = input$N0_max_test,
                       alpha0        = input$alpha0_test,
                       beta0         = input$beta0_test,
                       number_mcmc   = input$number_mcmc,
                       weibull_scale = input$weibull_scale_test,
                       weibull_shape = input$weibull_shape_test,
                       two_side      = input$two_side)
  })

  ### Fit control data
  posterior_control <- reactive({
    Binomial_posterior(y             = input$y_control,
                       N             = input$N_control,
                       y0            = input$y0_control,
                       N0            = input$N0_control,
                       N0_max        = input$N0_max_control,
                       alpha0        = input$alpha0_control,
                       beta0         = input$beta0_control,
                       number_mcmc   = input$number_mcmc,
                       weibull_scale = input$weibull_scale_control,
                       weibull_shape = input$weibull_shape_control,
                       two_side      = input$two_side)
  })

  ### Final function output
  f <- reactive({
    final(posterior_control(), posterior_test())
  })


  ### Output text
  output$text_test <- renderUI({

    if(posterior_test()$N0==0){
      out <- list(h3("Posterior test group information"),
                  wellPanel(p("No Prior Supplied")))
    } else{
      out <- list(h3("Posterior test group information"),
                  wellPanel(
                    p("Effective sample size of prior: ", posterior_test()$N0_effective),
                    p("Bayesian p-value (new vs historical data): ", posterior_test()$pvalue),
                    p("Loss function value: ", posterior_test()$alpha_loss),
                    p("N0 max: ", posterior_test()$N0_max)))
    }
    return(out)
  })

  output$text_control <- renderUI({
    if(posterior_control()$N0==0){
      out <- list(h3("Posterior control group information"),
                  wellPanel(p("No Prior Supplied")))
    } else{
      out <- list(h3("Posterior control group information"),
                  wellPanel(
                    p("Effective sample size of prior: ", posterior_control()$N0_effective),
                    p("Bayesian p-value (new vs historical data): ", posterior_control()$pvalue),
                    p("Loss function value: ", posterior_control()$alpha_loss),
                    p("N0 max: ", posterior_control()$N0_max)))
    }
    return(out)
  })

  ### Plot data for plots p and p1
  D <- reactive({
    D1 <- data.frame(information_sources = 'Posterior',
                     Group               = "Control",
                     y                   = f()$den_post_control$y,
                     x                   = f()$den_post_control$x)
    D2 <- data.frame(information_sources = "Current data",
                     Group               = "Control",
                     y                   = f()$den_flat_control$y,
                     x                   = f()$den_flat_control$x)
    D3 <- data.frame(information_sources = "Prior",
                     Group               = "Control",
                     y                   = f()$den_prior_control$y,
                     x                   = f()$den_prior_control$x)
    D4 <- data.frame(information_sources = 'Posterior',
                     Group               = "Test",
                     y                   = f()$den_post_test$y,
                     x                   = f()$den_post_test$x)
    D5 <- data.frame(information_sources = "Current data",
                     Group               = "Test",
                     y                   = f()$den_flat_test$y,
                     x                   = f()$den_flat_test$x)
    D6 <- data.frame(information_sources = "Prior",
                     Group               = "Test",
                     y                   = f()$den_prior_test$y,
                     x                   = f()$den_prior_test$x)

    if(input$N0_test==0 & input$N0_control==0){
      D <- as.data.frame(rbind(D1,D2,D4,D5))
    }
    if(input$N0_test==0 & input$N0_control!=0){
      D <- as.data.frame(rbind(D1,D2,D3,D4,D5))
    }
    if(input$N0_test!=0 & input$N0_control==0){
      D <- as.data.frame(rbind(D1,D2,D4,D5,D6))
    }
    if(input$N0_test!=0 & input$N0_control!=0){
      D <- as.data.frame(rbind(D1,D2,D3,D4,D5,D6,D6))
    }

    D$information_sources <- factor(D$information_sources,
                                    levels = (c("Posterior","Current data","Prior")))

    return(D)
  })


  ### P plot
  pPlot <- reactive({

    p <- ggplot(transform(D(), Group=factor(Group, levels=c("Test", "Control"))),
                aes(x=x,y=y)) +
      geom_line(size=1.5,aes(colour=information_sources,lty=information_sources)) +
      theme_bw() +
      facet_wrap(~Group, ncol=1, scales="free") +
      ylab("Density (PDF)") +
      xlab("Values")

    return(p)
  })

  output$plotP <- renderPlot({
    pPlot()
  })


  ### P1 plot
  p1Plot <- reactive({
    p1 <- ggplot(subset(D(),information_sources=="Posterior"),aes(x=x,y=y)) +
      geom_line(size=1.5,aes(colour=Group)) +
      ylab("Density (PDF)") +
      xlab("Values") +
      theme_bw()

    return(p1)
  })

  output$plotP1 <- renderPlot({
    p1Plot()
  })


  ### Plot histogram
  pHist <- reactive({
    ggplot(data.frame(x = f()$TestMinusControl_post),
           aes(x=x))  +
      geom_histogram(fill="dodgerblue", color="white") +
      xlab("Test-Control") +
      geom_vline(xintercept = 0, linetype = "longdash")
  })

  output$plotHist <- renderPlot({
    pHist()
  })


  ### Plot p3
  p3Plot <- reactive({

    if(input$two_side==1){
      p_value <- seq(0,1,,100)
      p_value <- ifelse(p_value>.5,1-p_value,p_value)
    }
    if(input$two_side==0){
      p_value <- seq(0,1,,100)
    }


    Loss_function_test <- pweibull(p_value,
                                   shape=posterior_test()$weibull_shape,
                                   scale=posterior_test()$weibull_scale)*posterior_test()$N0_max
    Loss_function_control <- pweibull(p_value,
                                      shape=posterior_control()$weibull_shape,
                                      scale=posterior_control()$weibull_scale)*posterior_control()$N0_max

    D1 <- data.frame(Group=c("Test"), y=Loss_function_test,x=seq(0,1,,100))
    D2 <- data.frame(Group=c("Test"), pvalue=c(posterior_test()$pvalue))
    D3 <- data.frame(Group=c("Test"), pvalue=c(posterior_test()$N0_effective))
    D4 <- data.frame(Group=c("Control"), y=Loss_function_control,x=seq(0,1,,100))
    D5 <- data.frame(Group=c("Control"), pvalue=c(posterior_control()$pvalue))
    D6 <- data.frame(Group=c("Control"), pvalue=c(posterior_control()$N0_effective))

    p3 <- ggplot()
    if(input$N0_test!=0){
      p3 <- p3 + geom_line(data=D1,aes(y=y,x=x,colour=Group),size=1) +
        geom_vline(data=D2, aes(xintercept=pvalue,colour=Group),lty=2) +
        geom_hline(data=D3, aes(yintercept=pvalue,colour=Group),lty=2)
    }

    if(input$N0_control!=0){
      p3 <- p3 + geom_line(data=D4,aes(y=y,x=x,colour=Group),size=1) +
        geom_vline(data=D5, aes(xintercept=pvalue,colour=Group),lty=2) +
        geom_hline(data=D6, aes(yintercept=pvalue,colour=Group),lty=2)
    }

    p3 <- p3 + facet_wrap(~Group, ncol=1) +
      theme_bw() +
      ylab("Effective sample size for historical data") +
      xlab("Bayesian p-value (new vs historical data)")

    return(p3)
  })

  output$plotP3 <- renderPlot({
    p3Plot()
  })

  ### Text for hypothesis test info
  output$hypothesis_text <- renderUI({
    out <- list(h3("Hypothesis testing information"),
                wellPanel(
                  p("We can define W as the difference between the event rates for the test group versus control group, i.e. W = test rate - control rate."),
                  p("Null Hypothesis (H_0): W ≥ 0.0%"),
                  p("Alternative Hypothesis (H_a): W < 0.0%"),
                  p("P(W < 0.0% | data) = ", mean(f()$TestMinusControl_post<0)),
                  p("We can accept H_a with a Probability of ", mean(f()$TestMinusControl_post<0))
                ))
    return(out)
  })


})
