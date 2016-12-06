################################################################################
# This code is used for estimating posterior samples from a binary outcome     #
# where an informative prior is used. The prior weight is determined using a   #
# loss function. In addition this code simulate many trials in order to get    #
# trial characteristics you must specify the parameters of the loss function   #
# as well as the maximum strength for the prior. This code assumes a           #
# non-adaptive trial.                                                          #
# This code is modeled after the methodologies developed by the MDIC working   #
# group: "Informing clinical trials using bench & simulations"                 #
# Section 1: of the code defines functions needed                              #
# Section 2: of the code estimates a posterior and the loss function value     #
#            given inputs                                                      #
# Section 3: of the code simulates the trial many times to get trial           #
#            characteristics                                                   #
# Developer: Tarek Haddad                                                      #
# Tarek.D.Haddad@Medtronic.com                                                 #
# Last modified:1/26/2016                                                      #
################################################################################


#############################################################
# This function produces the appropriate weight for your    #
# prior data assuming a binomial outcome                    #
# This function produces a scaler between 0 and 1           #
#############################################################
Loss_function <- function(y,N,y0,N0,alpha0,beta0,number_mcmc,weibull_shape,
                          weibull_scale, two_side){

  ### Theta for using Flat prior
  alpha_post_flat  <- y+alpha0
  beta_post_flat   <- N-y+beta0
  theta_post_flate <- rbeta(number_mcmc,alpha_post_flat,beta_post_flat) #flate prior

  ### Prior model
  alpha_prior <- y0 +alpha0
  beta_prior  <- N0-y0+beta0
  theta0      <- rbeta(number_mcmc, alpha_prior,beta_prior)#  flate prior

  ### Test of model vs real
  p_test <- mean(theta_post_flate<theta0)   # larger is higher failure

  ### Number of effective sample size given shape and scale loss function
  if (two_side == 0) {
    alpha_loss <- pweibull(p_test, shape = weibull_shape, scale = weibull_scale)
  } else if (two_side == 1){
    p_test1    <- ifelse(p_test > 0.5, 1 - p_test, p_test)
    alpha_loss <- pweibull(p_test1, shape = weibull_shape, scale = weibull_scale)
  }

  return(list(alpha_loss       = alpha_loss,
              pvalue           = p_test,
              theta_post_flate = theta_post_flate,
              theta0           = theta0))
}


#############################################################
# Calculates posterior estimation for Binomial distribution #
# Given alpha_loss value and maximum strength of(N0_max)    #
# prior if alpha_loss = 1                                   #
#############################################################
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


############################################################
##Combines the loss function and posterior              ####
###estimation into one function                         ####
############################################################
Binomial_posterior <- function(y,N,y0,N0,N0_max,alpha0,beta0,number_mcmc,
                               weibull_shape,weibull_scale,two_side){

  alpha_loss         <- Loss_function(y=y,N=N,y0=y0,N0=N0,alpha0=alpha0,beta0=beta0,number_mcmc=number_mcmc,
                                      weibull_shape=weibull_shape,weibull_scale=weibull_scale,
                                      two_side=two_side)
  Binomial_posterior <- theta_post_aug_bin(y=y,N=N,y0=y0,N0=N0,N0_max=N0_max,
                                           alpha_loss=alpha_loss$alpha_loss,
                                           alpha0=alpha0,
                                           beta0=beta0,number_mcmc=number_mcmc)

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


final <- function(posterior_control,posterior_test){
  den_post_control  <- density(posterior_control$Binomial_posterior,adjust = .5)
  den_flat_control  <- density(posterior_control$Binomial_posterior_flate,adjust = .5)
  den_prior_control <- density(posterior_control$Binomial_prior,adjust = .5)

  den_post_test     <- density(posterior_test$Binomial_posterior,adjust = .5)
  den_flat_test     <- density(posterior_test$Binomial_posterior_flate,adjust = .5)
  den_prior_test    <- density(posterior_test$Binomial_prior,adjust = .5)

  TestMinusControl_post <- posterior_test$Binomial_posterior-posterior_control$Binomial_posterior

  return(list(den_post_control      = den_post_control,
              den_flat_control      = den_flat_control,
              den_prior_control     = den_prior_control,
              den_post_test         = den_post_test,
              den_flat_test         = den_flat_test,
              den_prior_test        = den_prior_test,
              TestMinusControl_post = TestMinusControl_post))
}


