################################################################################
# This code is used for estimating posterior samples from a binary outcome     #
# where an informative prior is used. The prior weight is determined using a   #
# loss function. In addition this code simulate many trials in order to get    #
# trial characteristics you must specify the parameters of the loss function   #
# as well as the maximum strength for the prior. This code assumes a           #
# non-adaptive trial and is modeled after the methodologies developed by the   #
# MDIC working group: "Informing clinical trials using bench & simulations"    #
# Section 1: of the code defines functions needed                              #
# Section 2: of the code estimates a posterior and the loss function value     #
#            given inputs                                                      #
# Section 3: of the code simulates the trial many times to get trial           #
#            characteristics                                                   #
# Developer: Tarek Haddad                                                      #
# Tarek.D.Haddad@Medtronic.com                                                 #
# Last modified: 1/26/2016. Edited by Donnie Musgrove 7/14/2016                #
################################################################################

library(ggplot2)
library(parallel)

### This function produces appropriate weight for prior data assuming a
### binomial outcome. Produces a scalar between 0 and 1
Loss_function <- function(y,N,y0,N0,alpha0,beta0,number_mcmc,weibull_shape,
                          weibull_scale){

  ### Theta for using flat prior
  alpha_post_flat  <- y+alpha0
  beta_post_flat   <- N-y+beta0
  theta_post_flate <- rbeta(number_mcmc,alpha_post_flat,beta_post_flat) #flat prior

  #### Prior model
  alpha_prior <- y0 +alpha0
  beta_prior  <- N0-y0+beta0
  theta0      <- rbeta(number_mcmc,  alpha_prior,beta_prior) #flat prior

  ### Test of model vs real
  p_test <- mean(theta_post_flate<theta0)   #larger value == higher failure

  ### Number of effective sample size given shape and scale loss function
  alpha_loss <- pweibull(p_test,shape=weibull_shape,scale=weibull_scale)

  return(list(alpha_loss       = alpha_loss,
              pvalue           = p_test,
              theta_post_flate = theta_post_flate,
              theta0           = theta0))
}


### Calculate posterior estimation for Binomial distribution Given alpha_loss
### value and maximum strength of (N0_max) prior if alpha_loss==1
theta_post_aug_bin <- function(y,N,y0,N0,N0_max,alpha_loss,alpha0,beta0,
                               number_mcmc){

  effective_N0 <- N0_max*alpha_loss
  if(N0==0){
    alpha_prior <- alpha0
    beta_prior  <- beta0
  }else{
    alpha_prior <- (y0/N0)*effective_N0+alpha0
    beta_prior  <- effective_N0-(y0/N0)*effective_N0+beta0
  }
  alpha_post_aug <- y+alpha_prior
  beta_post_aug  <- N-y+beta_prior
  theta_post_aug <- rbeta(number_mcmc,alpha_post_aug,beta_post_aug)
  return(theta_post_aug)
}


### Combine loss function and posterior estimation into one function
Binomial_posterior <- function(y,N,y0,N0,N0_max,alpha0,beta0,number_mcmc,
                               weibull_shape,weibull_scale){

  alpha_loss <- Loss_function(y             = y,
                              N             = N,
                              y0            = y0,
                              N0            = N0,
                              alpha0        = alpha0,
                              beta0         = beta0,
                              number_mcmc   = number_mcmc,
                              weibull_shape = weibull_shape,
                              weibull_scale = weibull_scale)

  Binomial_posterior <- theta_post_aug_bin(y           = y,
                                           N           = N,
                                           y0          = y0,
                                           N0          = N0,
                                           N0_max      = N0_max,
                                           alpha_loss  = alpha_loss$alpha_loss,
                                           alpha0      = alpha0,
                                           beta0       = beta0,
                                           number_mcmc = number_mcmc)

  return(list(alpha_loss               = alpha_loss$alpha_loss,
              pvalue                   = alpha_loss$pvalue,
              Binomial_posterior       = Binomial_posterior,
              Binomial_posterior_flate = alpha_loss$theta_post_flate,
              Binomial_prior           = alpha_loss$theta0,
              weibull_scale            = weibull_scale,
              weibull_shape            = weibull_shape,
              y                        = y,
              N                        = N,
              y0                       = y0,
              N0                       = N0,
              N0_max                   = N0_max,
              N0_effective             = alpha_loss$alpha_loss*N0_max))
}


simulation_OPC <- function(confidence_level=0.95, true_rate, H0_rate, sim=100,
                           N, y0, N0, N0_max, alpha0=1, beta0=1,
                           number_mcmc=10000, weibull_scale=c(0.2,0.3),
                           weibull_shape=c(2,3)){

  runs      <- expand.grid(true_rate=true_rate,weibull_scale=weibull_scale,
                           weibull_shape=weibull_shape)

  trial_sim <- function(i,runs){
    y         <- rbinom(sim,N,runs$true_rate[i])
    post_rate <- rep(NA,sim)
    mean_N0   <- rep(NA,sim)
    for(j in 1:sim){
      posterior <- Binomial_posterior(y             = y[j],
                                      N             = N,
                                      y0            = y0,
                                      N0            = N0,
                                      N0_max        = N0_max,
                                      alpha0        = alpha0,
                                      beta0         = beta0,
                                      number_mcmc   = number_mcmc,
                                      weibull_scale = runs$weibull_scale[i],
                                      weibull_shape = runs$weibull_shape[i])

      post_rate[j] <- mean(posterior$Binomial_posterior<H0_rate)
      mean_N0[j]   <- posterior$N0_effective
    }

    Ha_pass_rate <- mean(post_rate>confidence_level)

    return(list(Ha_pass_rate = Ha_pass_rate,
                mean_N0      = mean(mean_N0)))
  }


  Ha_pass_rate1 <- mclapply(1:nrow(runs),FUN=trial_sim,runs=runs,mc.cores=1)
  Ha_pass_rate  <- rep(NA,length(Ha_pass_rate1))
  mean_N0       <- rep(NA,length(Ha_pass_rate1))

  for(jj in 1:length(Ha_pass_rate1)){
    Ha_pass_rate[jj] <- Ha_pass_rate1[[jj]]$Ha_pass_rate
    mean_N0[jj]      <- Ha_pass_rate1[[jj]]$mean_N0
  }

  avg_N0          <- round(mean_N0,1)
  historical_rate <- qbeta(.5,1+y0,1+N0-y0)
  results         <- data.frame(H0_rate=H0_rate, historical_rate, cbind(runs,avg_N0,Ha_pass_rate))

  return(results)
}

### Code execution
sim_results <- simulation_OPC(
  confidence_level = .95,
  true_rate = c(.01,.05,.075,.1),
  H0_rate = .1,
  sim = 50,
  N = 200,              #Number of  current subjects
  y0 = 15,              #Number of events observed  historical  data sets
  N0 = 200,             #Number of historical subjects
  N0_max = 200,         #Max effective sample size prior can receive when loss function equals 1
  alpha0 = 1,           #Noninformative Initial priors
  beta0 = 1,            #Noninformative Initial priors
  number_mcmc = 5000,   #Number of simulations to estimate posterior and loss function
  weibull_scale = c(.1,.2,.3), #Loss function parameter: location of a Weibull function
  weibull_shape = c(3,2))      #Loss function parameter: location of a Weibull function

sim_results$loss_function <- paste('scale=',sim_results$weibull_scale,",","shape=",sim_results$weibull_shape)

### Power curve
p3 <- ggplot(data=sim_results) +
  geom_line(aes(y=Ha_pass_rate,x=true_rate,colour=loss_function),size=1,lty=2) +
  geom_point(aes(y=Ha_pass_rate,x=true_rate,colour=loss_function),size=3) +
  geom_vline(aes(xintercept =historical_rate[1]),lty=2) +
  geom_text(aes(x=historical_rate[1], y=mean(Ha_pass_rate),label="Historical rate"))+
  ylab("Proportion of trials accepting HA") +
  xlab("True rate")


### Loss function plot
p_value <- seq(0,1,,100)
dd      <- subset(sim_results,true_rate==true_rate[1])
f <- function(i,dd,p_value){
  Loss_function <- pweibull(p_value,shape=dd$weibull_shape[i],scale=dd$weibull_scale[i])
  return(Loss_function)
}

ff <- sapply(1:nrow(dd),FUN=f,dd=dd,p_value=p_value)
pp3 <- ggplot() + ylab("Quantile") + xlab("Bayesian p-value")
for(i in 1:nrow(dd)){
  pp3 <- pp3 + geom_line(data=data.frame(p             = ff[,i],
                                         p_value       = p_value,
                                         loss_function = dd$loss_function[i]),
                         aes(y=p,x=p_value,colour=loss_function),
                         size=1,lty=2)
}

p3; pp3; sim_results
