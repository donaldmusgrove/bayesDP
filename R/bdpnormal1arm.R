#' bdpnormal1arm
#'
#' bdpnormal1arm
#'
#' @title bdpnormal1arm: bdpnormal1arm
#' @param mu_t numeric
#' @param sigma2_t numeric
#' @param N_t numeric
#' @param mu0_t numeric
#' @param sigma02_t numeric
#' @param N0_t numeric
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
#' @rdname bdpnormal1arm
#' @export bdpnormal1arm

################################################################################
# This code is used for estimating posterior samples from a Gaussian outcome   #
# where an informative prior is used. The prior weight is determined using a   #
# loss function. In addition this code simulate many trials in order to get    #
# trial characteristics you must specify the parameters of the loss function   #
# as well as the maximum strength for the prior. This code assumes a           #
# non-adaptive trial.                                                          #
# This code is modeled after the methodologies developed by the MDIC working   #
# group: "Informing clinical trials using bench & simulations"                 #
# Developer: Tarek Haddad                                                      #
# Tarek.D.Haddad@Medtronic.com                                                 #
# Last modified:1/26/2016                                                      #
################################################################################
#library(ggplot2)
#library(MCMCpack)
#library(survival)

#mu            = 10,
#sigma2        = 1,
#N             = 10,   #Number of  current subjects
#mu0           = 10,
#sigma02       = 1,
#N0            = 10,       #Number of historical subjects
#alpha_max     = 1,        #Max loss function weight
#number_mcmc   = 1000,     #Number of posterior simulations
#weibull_scale = 0.1,      #Loss function: Weibull location scale
#weibull_shape = 1,        #Loss function: Weibull location shape
#two_side      = two_side) # 0 == 1-sided, 1 === 2-sided

setGeneric("bdpnormal1arm",
           function(mu_t          = 10,
                    sigma2_t      = 1,
                    N_t           = 10,
                    mu0_t         = 10,
                    sigma02_t     = 1,
                    N0_t          = 10,
                    alpha_max     = 1,
                    weibull_scale = 0.1,
                    weibull_shape = 1,
                    number_mcmc   = 10000,
                    two_side      = 0,
                    inequality    = "<",
                    delta         = 0){
             standardGeneric("bdpnormal1arm")
           })

setMethod("bdpnormal1arm",
          signature(),
          function(mu_t          = 10,
                   sigma2_t      = 1,
                   N_t           = 10,
                   mu0_t         = 10,
                   sigma02_t     = 1,
                   N0_t          = 10,
                   alpha_max     = 1,
                   weibull_scale = 0.1,
                   weibull_shape = 1,
                   number_mcmc   = 10000,
                   two_side      = 0,
                   inequality    = "<",
                   delta         = 0){

################################################################################
### Functions
################################################################################

Loss_function <- function(mu, sigma2, N, mu0, sigma02, N0, alpha_max, number_mcmc,
                          weibull_shape, weibull_scale, two_side) {

  ### mu for using flat prior
  sigma2_post_flat <- rinvgamma(number_mcmc, (N-1)/2, ((N-1) * sigma2)/2)          #flat prior
  mu_post_flat     <- rnorm(number_mcmc, mu, (sigma2_post_flat/((N-1) + 1))^0.5)  #flat prior

  ### Prior model
  sigma2_post_flat0 <- rinvgamma(number_mcmc, (N0-1)/2, ((N0-1) * sigma02)/2)          #flat prior
  mu_post_flat0     <- rnorm(number_mcmc, mu0, (sigma2_post_flat0/((N0-1) + 1))^0.5)  #flat prior

  ### Test of model vs real
  p_test <- mean(mu_post_flat < mu_post_flat0)  # larger is higher failure

  ### Number of effective sample size given shape and scale loss function
  if (two_side == 0) {
    alpha_loss <- pweibull(p_test, shape = weibull_shape, scale = weibull_scale) * alpha_max
  } else if (two_side == 1) {
    p_test1    <- ifelse(p_test > 0.5, 1 - p_test, p_test)
    alpha_loss <- pweibull(p_test1, shape = weibull_shape, scale = weibull_scale) * alpha_max
  }

  return(list(alpha_loss    = alpha_loss,
              pvalue        = p_test,
              mu_post_flat  = mu_post_flat,
              mu0           = mu_post_flat0))
}


################################################################################
# Calculates posterior estimation for Binomial distribution given alpha_loss   #
################################################################################
mu_post_aug <- function(mu, sigma2, N, mu0, sigma02, N0, alpha_loss,
                        number_mcmc){
  effective_N0 <- N0 * alpha_loss
  sigma2_post  <- rinvgamma(number_mcmc, (N-1)/2, ((N-1) * sigma2)/2)    #flat prior
  sigma2_post0 <- rinvgamma(number_mcmc, (N0-1)/2, ((N0-1) * sigma02)/2) #flat prior

  mu1 <- (sigma2_post0 * N * mu + sigma2_post * effective_N0 * mu0)/(N *
                                                                       sigma2_post0 + sigma2_post * effective_N0)

  var_mu <- (sigma2_post * sigma2_post0)/(N * sigma2_post0 + sigma2_post *
                                            effective_N0)

  mu_post <- rnorm(number_mcmc, mu1, sqrt(var_mu))
  return(mu_post)
}


#################################################################################
### Combines the loss function and posterior estimation into one function
################################################################################
mu_posterior <- function(mu, sigma2, N, mu0, sigma02, N0, alpha_max, number_mcmc,
                         weibull_shape, weibull_scale, two_side) {

  alpha_loss <- Loss_function(mu            = mu,
                              sigma2        = sigma2,
                              N             = N,
                              mu0           = mu0,
                              sigma02       = sigma02,
                              N0            = N0,
                              alpha_max     = alpha_max,
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
                              alpha_loss  = alpha_loss$alpha_loss,
                              number_mcmc = number_mcmc)

  return(list(alpha_loss         = alpha_loss$alpha_loss,
              pvalue             = alpha_loss$pvalue,
              mu_posterior       = mu_posterior,
              mu_posterior_flat  = alpha_loss$mu_post_flat,
              mu_prior           = alpha_loss$mu0,
              weibull_scale      = weibull_scale,
              weibull_shape      = weibull_shape,
              mu                 = mu,
              N                  = N,
              mu0                = mu0,
              N0                 = N0,
              N0_effective       = alpha_loss$alpha_loss * N0))
}


final <- function(posterior_test) {
  den_post_test  <- density(posterior_test$mu_posterior, adjust = 0.5)
  den_flat_test  <- density(posterior_test$mu_posterior_flat, adjust = 0.5)
  den_prior_test <- density(posterior_test$mu_prior, adjust = 0.5)

  Testpost <- posterior_test$mu_posterior

  return(list(den_post_test  = den_post_test,
              den_flat_test  = den_flat_test,
              den_prior_test = den_prior_test,
              Testpost       = Testpost))
}


#results <- function(f,posterior_test,delta,two_side,inequality){





################################################################################
### Results
################################################################################

#two_side   <- 0    # 0 == 1-sided, 1 === 2-sided
#delta      <- 10   # delta value
#inequality <- ">"  # Inequality of alternate hypothesis

est <- mu_posterior(mu            = mu_t,
                    sigma2        = sigma2_t,
                    N             = N_t,
                    mu0           = mu0_t,
                    sigma02       = sigma02_t,
                    N0            = N0_t,
                    alpha_max     = alpha_max,
                    number_mcmc   = number_mcmc,
                    weibull_scale = weibull_scale,
                    weibull_shape = weibull_shape,
                    two_side      = two_side)

f1 <- final(posterior_test = est)

args1 <- list(mu_t          = mu_t,
              sigma2_t      = sigma2_t,
              N_t           = N_t,
              mu0_t         = mu0_t,
              sigma02_t     = sigma02_t,
              N0_t          = N0_t,
              alpha_max     = alpha_max,
              weibull_scale = weibull_scale,
              weibull_shape = weibull_shape,
              number_mcmc   = number_mcmc,
              two_side      = two_side,
              inequality    = inequality,
              delta         = delta)

#res1 <- results(f              = f1,
#                posterior_test = est,
#                delta          = delta,
#                two_side       = two_side,
#                inequality     = inequality)


### Plot outputs
#post_typeplot1 <- res1$post_typeplot

#densityplot1 <- res1$densityplot

#lossfun_plot1 <- res1$lossfun_plot

#lossfun_plot2 <- res1$lossfun_plot

### Text outputs
#hypothesis1 <- res1$hypothesis

#prior_for_test_group1 <- res1$prior_for_test_group

#setClass("bdpnormal1arm",
#         representation(post_typeplot1 = "ANY",
#                        densityplot1 = "ANY",
#                        lossfun_plot1 = "ANY",
#                        lossfun_plot2 = "ANY",
#                        hypothesis1 = "character",
#                        prior_for_test_group1 = "list"))

#me = new("bdpnormal1arm",
#         post_typeplot1 = post_typeplot1,
#         densityplot1 = densityplot1,
#         lossfun_plot1 = lossfun_plot1,
#         lossfun_plot2 = lossfun_plot2,
#         hypothesis1 = hypothesis1,
#         prior_for_test_group1 = prior_for_test_group1)

#me <- list(post_typeplot1 = post_typeplot1,
#           densityplot1 = densityplot1,
#           lossfun_plot1 = lossfun_plot1,
#           lossfun_plot2 = lossfun_plot2,
#           hypothesis1 = hypothesis1,
#           prior_for_test_group1 = prior_for_test_group1)

me <- list(est = est,
           f1 = f1,
           args1 = args1)

class(me) <- "bdpnormal1arm"

return(me)

})


#' plot
#'
#' plot
#'
#' @title plot: plot
#' @param x bdpnormal1arm
#'
#' @examples
#'
#' @rdname plot
#' @export plot
setMethod("plot", signature(x = "bdpnormal1arm"), function(x){

  f <- x$f1
  posterior_test <- x$est
  two_side <- x$args1$two_side
  inequality <- x$args1$inequality
  delta <- x$args1$delta

  D4 <- data.frame(information_sources = "Posterior",
                   group               = "Test",
                   y                   = f$den_post_test$y,
                   x                   = f$den_post_test$x)

  D5 <- data.frame(information_sources = "Current data",
                   group               = "Test",
                   y                   = f$den_flat_test$y,
                   x                   = f$den_flat_test$x)

  D6 <- data.frame(information_sources = "Prior",
                   group               = "Test",
                   y                   = f$den_prior_test$y,
                   x                   = f$den_prior_test$x)

  D <- as.data.frame(rbind(D4, D5, D6))

  D$information_sources <- factor(D$information_sources,
                                  levels = (c("Posterior", "Current data", "Prior")))

  post_typeplot <- ggplot(D, aes(x = x, y = y)) +
    geom_line(size = 2, aes(colour = information_sources, lty = information_sources)) +
    theme_bw() +
    facet_wrap(~group, ncol = 1, scale = "free") +
    ylab("Density (PDF)") +
    xlab("values")


  densityplot <- ggplot(subset(D, information_sources == "Posterior"), aes(x = x, y = y)) +
    geom_line(size = 2, aes(colour = group)) +
    ylab("Density (PDF)") +
    xlab("values") +
    theme_bw()

  if (two_side == 1) {
    p_value = seq(0, 1, , 100)
    p_value = ifelse(p_value > 0.5, 1 - p_value, p_value)
  }
  if (two_side == 0) {
    p_value = seq(0, 1, , 100)
  }

  Loss_function_test <- pweibull(p_value,
                                 shape = posterior_test$weibull_shape,
                                 scale = posterior_test$weibull_scale)*posterior_test$N0

  D1 <- data.frame(group = "test", y = Loss_function_test, x = seq(0, 1, , 100))
  D2 <- data.frame(group = c("test"), pvalue = c(posterior_test$pvalue))
  D3 <- data.frame(group = c("test"), pvalue = c(posterior_test$N0_effective))


  lossfun_plot <- ggplot() +
    geom_line(data = D1, aes(y = y, x = x, colour = group), size = 1) +
    geom_vline(data = D2, aes(xintercept = pvalue, colour = group), lty = 2) +
    geom_hline(data = D3, aes(yintercept = pvalue, colour = group), lty = 2) +
    facet_wrap(~group, ncol = 1) +
    theme_bw() +
    ylab("Effective sample size for historical data") +
    xlab("Bayesian p-value (new vs historical data)")

  op <- par(ask=TRUE)
  plot(post_typeplot)
  plot(densityplot)
  plot(lossfun_plot)
  par(op)
})

#' print
#'
#' print
#'
#' @title print: print
#' @param x bdpnormal1arm
#'
#' @examples
#'
#' @rdname print
#' @export print
setMethod("print", signature(x = "bdpnormal1arm"), function(x){

  f <- x$f1
  posterior_test <- x$est
  inequality <- x$args1$inequality
  delta <- x$args1$delta

  if (inequality == "<") {
    hypothesis <- paste("\"We can define mu as the mean for the test", "\n",
                        "Null Hypothesis (H_0): mu>", delta, "\n",
                        "Alternative Hypothesis (H_a): mu<", delta, "\n", "\n",
                        "P(mu<", delta, "|data)=", mean(f$Testpost < delta), "\n",
                        "We can accept H_a with a Probability of",
                        mean(f$Testpost < delta))
  }

  if (inequality == ">") {
    hypothesis <- paste("\"Define mu as the mean of the data", "\n",
                        "Null Hypothesis (H_0): mu<", delta, "\n",
                        "Alternative Hypothesis (H_a): mu>", delta, "\n", "\n",
                        "P(mu>", delta, "|data)=", mean(f$Testpost > delta), "\n",
                        "We can accept H_a with a Probability of",
                        mean(f$Testpost > delta))
  }

  ### Print
  prior_for_test_group <- list(`Effective sample size of prior (for test group)` = posterior_test$N0_effective,
                               `Bayesian p-value (new vs historical data)`       = posterior_test$pvalue,
                               `Loss function value`                             = posterior_test$alpha_loss,
                               `Sample size of prior (for test group)`           = posterior_test$N0)

  print(cat(hypothesis))
  print(prior_for_test_group)
})

#' summary
#'
#' summary
#'
#' @title summary: summary
#' @param object bdpnormal1arm
#'
#' @examples
#'
#' @rdname summary
#' @export summary
setMethod("summary", signature(object = "bdpnormal1arm"), function(object){

  f <- object$f1
  posterior_test <- object$est
  inequality <- object$args1$inequality
  delta <- object$args1$delta

  if (inequality == "<") {
    hypothesis <- paste("\"We can define mu as the mean for the test", "\n",
                        "Null Hypothesis (H_0): mu>", delta, "\n",
                        "Alternative Hypothesis (H_a): mu<", delta, "\n", "\n",
                        "P(mu<", delta, "|data)=", mean(f$Testpost < delta), "\n",
                        "We can accept H_a with a Probability of",
                        mean(f$Testpost < delta))
  }

  if (inequality == ">") {
    hypothesis <- paste("\"Define mu as the mean of the data", "\n",
                        "Null Hypothesis (H_0): mu<", delta, "\n",
                        "Alternative Hypothesis (H_a): mu>", delta, "\n", "\n",
                        "P(mu>", delta, "|data)=", mean(f$Testpost > delta), "\n",
                        "We can accept H_a with a Probability of",
                        mean(f$Testpost > delta))
  }

  ### summary
  prior_for_test_group <- list(`Effective sample size of prior (for test group)` = posterior_test$N0_effective,
                               `Bayesian p-value (new vs historical data)`       = posterior_test$pvalue,
                               `Loss function value`                             = posterior_test$alpha_loss,
                               `Sample size of prior (for test group)`           = posterior_test$N0)

  print(cat(hypothesis))
  print(prior_for_test_group)
  argsdf <- data.frame(t(data.frame(object$args1)))
  names(argsdf) <- "args"
  print(argsdf)
})
