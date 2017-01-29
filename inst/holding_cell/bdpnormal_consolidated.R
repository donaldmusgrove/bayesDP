#' This code is used for estimating posterior samples from a Gaussian outcome
#' where an informative prior is used. The prior weight is determined using a
#' loss function. In addition this code simulate many trials in order to get
#' trial characteristics you must specify the parameters of the loss function
#' as well as the maximum strength for the prior. This code assumes a
#' non-adaptive trial.
#' This code is modeled after the methodologies developed by the MDIC working
#' group: "Informing clinical trials using bench & simulations"
#' Developer: Tarek Haddad
#' Tarek.D.Haddad@Medtronic.com
#' Last modified:1/26/2016
#'
#'
#' bdpnormal2arm
#'
#' @title bdpnormal2arm: bdpnormal2arm
#' @param mu_t numeric
#' @param sigma_t numeric
#' @param N_t numeric
#' @param mu_c numeric
#' @param sigma_c numeric
#' @param N_c numeric
#' @param mu0_t numeric
#' @param sigma0_t numeric
#' @param N0_t numeric
#' @param mu0_c numeric
#' @param sigma0_c numeric
#' @param N0_c numeric
#' @param alpha_max numeric
#' @param weibull_scale numeric
#' @param weibull_shape numeric
#' @param number_mcmc numeric
#' @param two_side numeric
#' @param inequality character
#' @param delta numeric
#'
#' @examples
#'
#' @rdname bdpnormal2arm
#' @export bdpnormal2arm

setGeneric("bdpnormal2arm",
           function(mu_t = NULL,
                    sigma_t = NULL,
                    N_t = NULL,
                    mu_c = NULL,
                    sigma_c = NULL,
                    N_c = NULL,
                    mu0_t = NULL,
                    sigma0_t = NULL,
                    N0_t = NULL,
                    mu0_c = NULL,
                    sigma0_c = NULL,
                    N0_c = NULL,  # up to here null
                    alpha_max = 1, # default 1
                    weibull_scale = 0.135, #  0.135
                    weibull_shape = 3, # 3
                    number_mcmc  = 10000, # 10000 good
                    two_side = 1){ # get rid of this
             standardGeneric("bdpnormal2arm")
           })

setMethod("bdpnormal2arm",
          signature(),
          function(mu_t = NULL,
                   sigma_t = NULL,
                   N_t = NULL,
                   mu_c = NULL,
                   sigma_c = NULL,
                   N_c = NULL,
                   mu0_t = NULL,
                   sigma0_t = NULL,
                   N0_t = NULL,
                   mu0_c = NULL,
                   sigma0_c = NULL,
                   N0_c = NULL,  # up to here null
                   alpha_max = 1, # default 1
                   weibull_scale = 0.135, #  0.135
                   weibull_shape = 3, # 3
                   number_mcmc  = 10000, # 10000 good
                   two_side = 1){ # get rid of this

################################################################################
# Produce prior data weight (scalar between 0 and 1) assuming a mu outcome     #
################################################################################
Loss_function <- function(mu, sigma, N, mu0, sigma0, N0, alpha_max, number_mcmc,
                          weibull_shape, weibull_scale, two_side){

  ### mu for using flat prior
  sigma2_post_flat <- rinvgamma(number_mcmc, (N - 1)/2, ((N - 1) * sigma^2)/2)
  mu_post_flat     <- rnorm(number_mcmc, mu, (sigma2_post_flat/((N-1) + 1))^0.5)

  ### Prior model (flat priors)
  sigma2_post_flat0 <- rinvgamma(number_mcmc, (N0-1)/2, ((N0-1) * sigma0^2)/2)
  mu_post_flat0     <- rnorm(number_mcmc, mu0, (sigma2_post_flat0/((N0-1) + 1))^0.5)

  ### Test of model vs real
  p_test <- mean(mu_post_flat < mu_post_flat0)  # larger is higher failure

  ### Number of effective sample size given shape and scale loss function
  if (two_side == 0) {
    alpha_loss <- pweibull(p_test, shape = weibull_shape, scale = weibull_scale) * alpha_max
  } else if (two_side == 1){
    p_test1    <- ifelse(p_test > 0.5, 1 - p_test, p_test)
    alpha_loss <- pweibull(p_test1, shape = weibull_shape, scale = weibull_scale) * alpha_max
  }

  return(list(alpha_loss   = alpha_loss,
              pvalue       = p_test,
              mu_post_flat = mu_post_flat,
              mu0          = mu_post_flat0))
}
