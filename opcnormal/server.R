library(ggplot2)
library(knitr)
library(MCMCpack)
library(parallel)
library(rmarkdown)
library(shiny)
library(shinyBS)
library(survival)

shinyServer(function(input, output) {
  source("post_estimation_normal_functions.R")
  source("Normal_OPC_shiny.R")

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

