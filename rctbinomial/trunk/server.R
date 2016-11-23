


### Name: Bayesian loss function RCT

    
shinyServer(function(input, output, session){

  ### Load external R functions based on distribution choice
  #observeEvent(input$distribution,{
  #  if(input$distribution %in% c("binomial","a1","a2","a3")){
    source("www/post_estimation_binomial_functions.R")
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










