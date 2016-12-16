################################################################################
# This code is used for estimating posterior samples from a Gaussian outcome   #
# where an informative prior is used. The prior weight is determined using a   #
# loss function. In addition this code simulate many trials in order to get    #
# trial characteristics you must specify the parameters of the loss function   #
# as well as the maximum strength for the prior. This code assumes a           #
# non-adaptive trial. This code is modeled after the methodologies developed   #
# by the MDIC working group:                                                   #
# 'Informing clinical trials using bench & simulations'                        #
# Developer: Tarek Haddad                                                      #
# Tarek.D.Haddad@Medtronic.com                                                 #
# Last modified: 9/10/2016                                                     #
################################################################################


################################################################################
### Results
################################################################################

## declare the display generic
setGeneric("rtcnormal", function(f, ...)
  standardGeneric("rtcnormal")
)

setMethod("rtcnormal",
          signature(f = "ANY"),
          function(f){
            message("Wrong object")
          })

setMethod("rtcnormal",
          signature(f = "missing"),
          function(f){
            message("Missing object")
          })

setMethod("rtcnormal",
          signature(f = "numeric"),
results <- function(f,posterior_test,posterior_control,two_side,inequality,
                    N0_t,N0_c,delta=2){
  D1 <- data.frame(information_sources='Posterior',
                   group="Control",
                   y=f$den_post_control$y,
                   x=f$den_post_control$x)
  D2 <- data.frame(information_sources="Current data",
                   group="Control",
                   y=f$den_flat_control$y,
                   x=f$den_flat_control$x)
  D3 <- data.frame(information_sources="Prior",
                   group="Control",
                   y=f$den_prior_control$y,
                   x=f$den_prior_control$x)
  D4 <- data.frame(information_sources='Posterior',
                   group="Test",
                   y=f$den_post_test$y,
                   x=f$den_post_test$x)
  D5 <- data.frame(information_sources="Current data",
                   group="Test",
                   y=f$den_flat_test$y,
                   x=f$den_flat_test$x)
  D6 <- data.frame(information_sources="Prior",
                   group="Test",
                   y=f$den_prior_test$y,
                   x=f$den_prior_test$x)

  if(N0_t==0 & N0_c==0){
    D <- as.data.frame(rbind(D4,D5,D1,D2))
  }
  if(N0_t==0 & N0_c!=0){
    D <- as.data.frame(rbind(D4,D5,D1,D2,D3))
  }
  if(N0_t!=0 & N0_c==0){
    D <- as.data.frame(rbind(D4,D5,D6,D1,D2))
  }
  if(N0_t!=0 & N0_c!=0){
    D <- as.data.frame(rbind(D4,D5,D6,D1,D2,D3))
  }

  D$information_sources <- factor(D$information_sources,
                                  levels = (c("Posterior","Current data","Prior")))

  post_typeplot <- ggplot(D,aes(x=x,y=y)) +
    geom_line(size=2,aes(colour=information_sources,lty=information_sources)) +
    theme_bw() +
    facet_wrap(~group, ncol=1,scale='free') +
    ylab("Density (PDF)") +
    xlab("values")

  densityplot <- ggplot(subset(D,information_sources=="Posterior"),
                        aes(x=x,y=y)) +
    geom_line(size=2,aes(colour=group)) +
    ylab("Density (PDF)") +
    xlab("values") +
    theme_bw()


  if(two_side==1){
    p_value <- seq(0,1,,100)
    p_value <- ifelse(p_value>.5,1-p_value,p_value)
  }
  if(two_side==0){
    p_value <- seq(0,1,,100)
  }

  Loss_function_test <- pweibull(p_value,
                                 shape=posterior_test$weibull_shape,
                                 scale=posterior_test$weibull_scale)*posterior_test$N0_max

  Loss_function_control <- pweibull(p_value,
                                    shape=posterior_control$weibull_shape,
                                    scale=posterior_control$weibull_scale)*posterior_control$N0_max

  D1 <- data.frame(group="test",y=Loss_function_test,x=seq(0,1,,100))
  D2 <- data.frame(group=c("test"),pvalue=c(posterior_test$pvalue))
  D3 <- data.frame(group=c("test"),pvalue=c(posterior_test$N0_effective))


  D4 <- data.frame(group="control",y=Loss_function_control,x=seq(0,1,,100))
  D5 <- data.frame(group=c("control"),pvalue=c(posterior_control$pvalue))
  D6 <- data.frame(group=c("control"),pvalue=c(posterior_control$N0_effective))


  lossfun_plot <- ggplot()
  if(N0_t!=0){
    lossfun_plot <- lossfun_plot +
      geom_line(data=D1,aes(y=y,x=x,colour=group),size=1) +
      geom_vline(data=D2, aes(xintercept =pvalue,colour=group),lty=2) +
      geom_hline(data=D3, aes(yintercept =pvalue,colour=group),lty=2)
  }
  if(N0_c!=0){
    lossfun_plot  <- lossfun_plot +
      geom_line(data=D4,aes(y=y,x=x,colour=group),size=1) +
      geom_vline(data=D5, aes(xintercept =pvalue,colour=group),lty=2) +
      geom_hline(data=D6, aes(yintercept =pvalue,colour=group),lty=2)
  }

  lossfun_plot <- lossfun_plot +
    facet_wrap(~group, ncol=1) +
    theme_bw() +
    ylab("Effective sample size for historical data") +
    xlab("Bayesian p-value (new vs historical data)")

  if(inequality=="<"){
  hypothesis <- paste('"We can define W as the difference between the means
                      for the test group versus control group, i.e. W=test
                      mean - control mean', '\n', 'Null Hypothesis (H_0): W>',
                      delta, '\n', "Alternative Hypothesis (H_a): W<", delta,
                      '\n', '\n', "P(W<",delta,"|data)=",
                      mean(f$TestMinusControl_post<delta),
                      '\n', "We can accept H_a with a Probability of",
                      mean(f$TestMinusControl_post<delta))
  } else if(inequality==">"){
  hypothesis <- paste('"We can define W as the difference between the means
                      for the test group versus control group, i.e. W=test
                      mean - control mean', '\n', 'Null Hypothesis (H_0): W<',
                      delta, '\n', "Alternative Hypothesis (H_a): W>", delta,
                      '\n', '\n', "P(W<",delta,"|data)=",
                      mean(f$TestMinusControl_post>delta),
                      '\n', "We can accept H_a with a Probability of",
                      mean(f$TestMinusControl_post>delta))
  }

  if(posterior_test$N0==0){
    prior_for_test_group <- "No Prior Supplied"
  } else{
    prior_for_test_group <- list("Effective sample size of prior(for test group)"=posterior_test$N0_effective,
                                 "Bayesian p-value (new vs historical data)"=posterior_test$pvalue,
                                 'loss function value'=posterior_test$alpha_loss,
                                 'N0_max'=posterior_test$N0_max)
  }

  if(posterior_control$N0==0){
    prior_for_control_group <- "No Prior Supplied"
  } else{
    prior_for_control_group <- list("Effective sample size of prior(for control group)"=posterior_control$N0_effective,
                                    "Bayesian p-value (new vs historical data)"=posterior_control$pvalue,
                                    'loss function value'=posterior_control$alpha_loss,
                                    'N0_max'=posterior_control$N0_max)
  }

  return(list(prior_for_test_group    = prior_for_test_group,
              prior_for_control_group = prior_for_control_group,
              post_typeplot           = post_typeplot,
              densityplot             = densityplot,
              lossfun_plot            = lossfun_plot,
              hypothesis              = hypothesis))
}
)
