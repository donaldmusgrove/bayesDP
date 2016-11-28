
```{r}

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

library(ggplot2)
library(MCMCpack)
library(survival)


################################################################################ 
### Functions 
################################################################################
### Produce prior data weight (scalar between 0 and 1) assuming a mu outcome
Loss_function <- function(mu, sigma2, N, mu0, sigma02, N0, number_mcmc, 
                          weibull_shape, weibull_scale, two_side){

  ### mu for using flat prior
  sigma2_post_flate <- rinvgamma(number_mcmc, (N - 1)/2, ((N - 1) * sigma2)/2)
  mu_post_flate     <- rnorm(number_mcmc, mu, (sigma2_post_flate/((N-1)+1))^0.5) 
  
  ### Prior model (flat priors)
  sigma2_post_flate0 <- rinvgamma(number_mcmc, (N0-1)/2, ((N0-1)*sigma02)/2) 
  mu_post_flate0     <- rnorm(number_mcmc, mu0, 
                              (sigma2_post_flate0/((N0-1)+1))^0.5) 
  
  ### Test of model vs real
  p_test <- mean(mu_post_flate < mu_post_flate0)  # larger is higher failure
  
  ### Number of effective sample size given shape and scale loss function
  if (two_side == 0) {
    alpha_loss <- pweibull(p_test, shape = weibull_shape, scale = weibull_scale)
  } else if (two_side == 1){
    p_test1    <- ifelse(p_test > 0.5, 1 - p_test, p_test)
    alpha_loss <- pweibull(p_test1, shape = weibull_shape, scale = weibull_scale)
  }
  return(list(alpha_loss    = alpha_loss, 
              pvalue        = p_test, 
              mu_post_flate = mu_post_flate, 
              mu0           = mu_post_flate0))
}


### Estimate posterior for mu given alpha_loss value and maximum strength of(N0_max) 
### Prior if alpha_loss=1
mu_post_aug <- function(mu, sigma2, N, mu0, sigma02, N0, N0_max, alpha_loss, 
                        number_mcmc) {
  if (N0 != 0){
    effective_N0 <- N0_max * alpha_loss
    sigma2_post  <- rinvgamma(number_mcmc, (N-1)/2, ((N-1)*sigma2)/2)
    sigma2_post0 <- rinvgamma(number_mcmc, (N0-1)/2, ((N0-1)*sigma02)/2) 
    
    mu1 <- (sigma2_post0*N*mu + sigma2_post*effective_N0*mu0)/(N*sigma2_post0 + 
            sigma2_post*effective_N0)
    var_mu <- (sigma2_post*sigma2_post0)/(N*sigma2_post0 + 
               sigma2_post*effective_N0)
  } else {
    var_mu <- rinvgamma(number_mcmc, (N - 1)/2, ((N - 1) * sigma2)/2)
    mu1    <- rnorm(number_mcmc, mu, (var_mu/((N - 1) + 1))^0.5) 
    
  }
  mu_post <- rnorm(number_mcmc, mu1, sqrt(var_mu))
  return(mu_post)
}


### Combine loss function and posterior estimation into one function
mu_posterior <- function(mu, sigma2, N, mu0, sigma02, N0, N0_max, number_mcmc, 
                         weibull_shape, weibull_scale, two_side) {
  if (N0 != 0) {
    alpha_loss <- Loss_function(mu            = mu, 
                                sigma2        = sigma2, 
                                N             = N, 
                                mu0           = mu0, 
                                sigma02       = sigma02, 
                                N0            = N0, 
                                number_mcmc   = number_mcmc, 
                                weibull_shape = weibull_shape, 
                                weibull_scale = weibull_scale, 
                                two_side      = two_side)
    
    mu_posterior <- mu_post_aug(mu          = mu, 
                                sigma2      = sigma2, 
                                N           = N, 
                                mu0         = mu0, 
                                sigma02     = sigma02, 
                                N0          = N0, 
                                N0_max      = N0_max, 
                                alpha_loss  = alpha_loss$alpha_loss, 
                                number_mcmc = number_mcmc)
  } else {
    mu_posterior <- mu_post_aug(mu          = mu, 
                                sigma2      = sigma2, 
                                N           = N, 
                                mu0         = mu0, 
                                sigma02     = sigma02, 
                                N0          = N0, 
                                N0_max      = N0_max, 
                                alpha_loss  = 0, 
                                number_mcmc = number_mcmc)
      
    alpha_loss <- list(alpha_loss    = 0, 
                       pvalue        = 0, 
                       mu0           = rnorm(100), 
                       mu_post_flate = mu_posterior)
  }
  return(list(alpha_loss         = alpha_loss$alpha_loss, 
              pvalue             = alpha_loss$pvalue, 
              mu_posterior       = mu_posterior, 
              mu_posterior_flate = alpha_loss$mu_post_flate, 
              mu_prior           = alpha_loss$mu0, 
              weibull_scale      = weibull_scale, 
              weibull_shape      = weibull_shape, 
              mu                 = mu, 
              N                  = N, 
              mu0                = mu0, 
              N0                 = N0, 
              N0_max             = N0_max, 
              N0_effective       = 0))
}

final <- function(posterior_control, posterior_test) {
  den_post_control  <- density(posterior_control$mu_posterior, adjust = 0.5)
  den_flat_control  <- density(posterior_control$mu_posterior_flate, 
                               adjust = 0.5)
  den_prior_control <- density(posterior_control$mu_prior, adjust = 0.5)
  
  den_post_test  <- density(posterior_test$mu_posterior, adjust = 0.5)
  den_flat_test  <- density(posterior_test$mu_posterior_flate, adjust = 0.5)
  den_prior_test <- density(posterior_test$mu_prior, adjust = 0.5)
  
  TestMinusControl_post <- posterior_test$mu_posterior - posterior_control$mu_posterior
  
  return(list(den_post_control      = den_post_control, 
              den_flat_control      = den_flat_control, 
              den_prior_control     = den_prior_control, 
              den_post_test         = den_post_test, 
              den_flat_test         = den_flat_test, 
              den_prior_test        = den_prior_test, 
              TestMinusControl_post = TestMinusControl_post))
}

final1 <- function(posterior_test) {
  den_post_test  <- density(posterior_test$mu_posterior, adjust = 0.5)
  den_flat_test  <- density(posterior_test$mu_posterior_flate, adjust = 0.5)
  den_prior_test <- density(posterior_test$mu_prior, adjust = 0.5)
  
  Testpost <- posterior_test$mu_posterior
  
  return(list(den_post_test  = den_post_test, 
              den_flat_test  = den_flat_test, 
              den_prior_test = den_prior_test, 
              Testpost       = Testpost))
}


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


################################################################################ 
### Results 
################################################################################

two_side   <- 1   # 0 == 1-sided, 1 === 2-sided
inequality <- "<" # Inequality of alternate hypothesis
N0_t       <- 10  #Number of historical subjects in test group
N0_c       <- 0   #Number of historical subjects in control group 
delta      <- 0   #Non-inferiority zone value

posterior_test <- mu_posterior(
  mu            = 10,    #Number of events observed  current data sets   
  sigma2        = 5,
  N             = 10,    #Number of  current subjects                        
  mu0           = 9, 
  sigma02       = 5 ,    #Number of events observed  historical  data sets   
  N0            = N0_t,  #Number of historical subjects       
  N0_max        = 20,    #Maximum effective sample size prior can receive when the loss function equals 1
  number_mcmc   = 10000, #Number of simulations to estimate posterior and loss function
  weibull_scale = .2,    #Loss function parameter controlling the location of a weibull function
  weibull_shape = 2,     #Loss function parameter controlling the location of a weibull function
  two_side      = two_side)

    
posterior_control <- mu_posterior(
  mu            = 10,    #Number of events observed  current data sets   
  sigma2        = 2,
  N             = 20,    #Number of  current subjects                        
  mu0           = 12, 
  sigma02       = 1,     #Number of events observed  historical  data sets   
  N0            = N0_c,  #Number of historical subjects       
  N0_max        = 20,    #The maximum effective sample size the prior can receive when the loss function equals 1
  number_mcmc   = 10000, #Number of simulations to estimate posterior and loss function
  weibull_scale = .2,    #Loss function parameter controlling the location of a weibull function
  weibull_shape = 2,     #Loss function parameter controlling the location of a weibull function
  two_side      = two_side)

  
  
f1 <- final(posterior_test    = posterior_test,
            posterior_control = posterior_control)

  
res1 <- results(f                 = f1,
                posterior_test    = posterior_test,
                posterior_control = posterior_control,
                two_side          = two_side,
                inequality        = inequality,
                N0_t              = N0_t,
                N0_c              = N0_c,
                delta             = delta)


### Plot outputs
post_typeplot1 <- res1$post_typeplot
densityplot1   <- res1$densityplot
lossfun_plot1  <- res1$lossfun_plot


### Text outputs
hypothesis1              <- res1$hypothesis
prior_for_test_group1    <- res1$prior_for_test_group
prior_for_control_group1 <- res1$prior_for_control_group


### Display outputs
post_typeplot1
densityplot1 
lossfun_plot1
cat(hypothesis1)
prior_for_test_group1 
prior_for_control_group1

```
