

    
shinyServer(function(input, output, session){

  ### Set the number of cores, conditional on platform
  mcCores <- ifelse(.Platform$OS.type=="windows",1,20)


  ### Load external R functions based on distribution choice
  observeEvent(input$distribution,{
    if(input$distribution %in% c("binomial","a1","a2","a3")){
    source("www/binomial.R")
    }
  })
  
  
  ### Parse vector inputs into useable form
  inputVec <- eventReactive(input$button, {
    true_rate     <- as.numeric(unlist(strsplit(input$true_rate,",")))
    weibull_scale <- as.numeric(unlist(strsplit(input$weibull_scale,",")))
    weibull_shape <- as.numeric(unlist(strsplit(input$weibull_shape,",")))
    
    list(true_rate     = true_rate,
         weibull_scale = weibull_scale,
         weibull_shape = weibull_shape)
  })
   
  
  sim_results <- eventReactive(input$button, {
  
    withProgress(message="===> Progress: Carrying out simulations",value=0.6,{
  
    res <- simulation_OPC(confidence_level = input$CL, 
                          true_rate        = inputVec()$true_rate,
                          H0_rate          = input$H0_rate, 
                          sim              = input$sim, 
                          N                = input$N, 
                          y0               = input$y0,
                          N0               = input$N0,
                          N0_max           = input$N0_max,
                          alpha0           = input$alpha0, 
                          beta0            = input$beta0, 
                          number_mcmc      = input$number_mcmc,
                          weibull_scale    = inputVec()$weibull_scale, 
                          weibull_shape    = inputVec()$weibull_shape,
                          mcCores          = mcCores)
                          
    res$loss_function <- paste0('scale=', res$weibull_scale,", ","shape=",res$weibull_shape)
    
    setProgress(1)
    })
    
    res
  })
  
  
  p3 <- reactive({
    ggplot(data=sim_results()) +
      geom_line(aes(y=Ha_pass_rate, x=true_rate, colour=loss_function), size=1,lty=2) +
      geom_point(aes(y=Ha_pass_rate, x=true_rate, colour=loss_function), size=3) +
      geom_vline(aes(xintercept=historical_rate[1]),lty=2)+
      geom_text(aes(x=historical_rate[1], y=mean(Ha_pass_rate), label="Historical rate")) +
      ylab("Proportion of trials accepting HA") +
      xlab("True rate")
  })
  

  pp3 <- reactive({
    p_value <- seq(0,1,,100)
    dd      <- subset(sim_results(),true_rate==true_rate[1])

    f <- function(i,dd,p_value){
      Loss_function <- pweibull(p_value, shape=dd$weibull_shape[i], scale=dd$weibull_scale[i])
      return(Loss_function)
    }

    ff  <- sapply(1:nrow(dd),FUN=f,dd=dd,p_value=p_value)
    pp3 <- ggplot() + ylab("Quantile") + xlab("Bayesian p-value")
    for(i in 1:nrow(dd)){
      pp3 <- pp3 + geom_line(data=data.frame(p             = ff[,i], 
                                             p_value       = p_value, 
                                             loss_function = dd$loss_function[i]),
                             aes(y = p, x = p_value, colour = loss_function),
                             size=1, lty=2)
    } 
    pp3
  })
  
  
  sim_res <- reactive({
    if(is.null(sim_results())){
      return(NULL)
    } else{
      res <- sim_results()
      names(res) <- c('H0 Rate', 'Historical Rate', 'True Rate', 'Weibull Scale', 'Weibull Shape',	 
                      'Avg N0', 'HA Pass Rate', 'Loss Function')
                      
      res$`Historical Rate` <- round(res$`Historical Rate`, 3)
      return(res[,-c(4,5)])
    }
  
  })
  
  
  output$sim_results <- renderTable({
    sim_res()
  })
  
  output$p3 <- renderPlot({
    p3()
  })

  output$pp3 <- renderPlot({
    pp3()
  })
   
  
  ### Download a report  
  output$downloadReport <- downloadHandler(
    filename = function() {
      paste('BayesianLossFunctionReport', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'
      ))
    },  
    
    content = function(file) {
      src <- normalizePath('report.Rmd')

      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'report.Rmd')
          
      out <- render('report.Rmd', 
                    switch(input$format,
                           PDF  = pdf_document(), 
                           HTML = html_document(), 
                           Word = word_document()
                           ))     
      file.rename(out, file)
    }
  )
  
  
})