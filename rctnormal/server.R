library(shiny)
library(survival)
library(ggplot2)
library(MCMCpack)
library(knitr)
library(rmarkdown)


shinyServer(function(input, output) {
  source("post_estimation_normal_functions.R")
  source("Normal_RCT_shiny.R")


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

