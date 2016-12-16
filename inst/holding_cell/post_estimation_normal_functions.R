
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



