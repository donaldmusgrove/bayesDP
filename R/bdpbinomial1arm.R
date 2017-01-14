#' bdpbinomial1arm
#'
#' bdpbinomial1arm
#'
#' @title bdpbinomial1arm: bdpbinomial1arm
#' @param y numeric
#' @param N numeric
#' @param y0 numeric
#' @param N0 numeric
#' @param alpha_max numeric
#' @param a0 numeric
#' @param b0 numeric
#' @param number_mcmc numeric
#' @param weibull_shape numeric
#' @param weibull_scale numeric
#' @param H0 numeric
#' @param two_side character
#' @param inequality character
#'
#' @examples
#'
#' @rdname bdpbinomial1arm
#' @export bdpbinomial1arm


################################################################################
# This code is used for estimating posterior samples from a binary outcome     #
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

setGeneric("bdpbinomial1arm",
           function(y             = 1,     #n events: current
                    N             = 400,   #n subjects: current
                    y0            = 10,    #n events: historical
                    N0            = 100,   #n subjects: historical
                    alpha_max     = 1,     #Max loss function weight
                    a0            = 1,     #Noninformative Initial priors
                    b0            = 1,     #Noninformative Initial priors
                    number_mcmc   = 10000, #n simulations
                    weibull_scale = .05,   #Loss function scale parameter
                    weibull_shape = 2,     #Loss function shape
                    H0 = 10,
                    two_side      = 0,
                    inequality = "<"){
             standardGeneric("bdpbinomial1arm")
           })

setMethod("bdpbinomial1arm",
          signature(y = "numeric"),
          function(y             = 1,     #n events: current
                   N             = 400,   #n subjects: current
                   y0            = 10,    #n events: historical
                   N0            = 100,   #n subjects: historical
                   alpha_max     = 1,     #Max loss function weight
                   a0            = 1,     #Noninformative Initial priors
                   b0            = 1,     #Noninformative Initial priors
                   number_mcmc   = 10000, #n simulations
                   weibull_scale = .05,   #Loss function scale parameter
                   weibull_shape = 2,     #Loss function shape parameter
                   H0 = 10,
                   two_side      = 0,
                   inequality = "<"){


################################################################################
# Estimate weight for prior data assuming a binomial outcome                   #
################################################################################
Loss_function <- function(y, N, y0, N0, alpha_max, a0, b0, number_mcmc, weibull_shape,
                          weibull_scale, two_side=0){

  ### Theta for using flat prior
  a_post_flat     <- y + a0
  b_post_flat     <- N - y + b0
  theta_post_flat <- rbeta(number_mcmc, a_post_flat, b_post_flat)

  ### Prior model
  a_prior  <- y0 + a0
  b_prior  <- N0 - y0 + b0
  theta0   <- rbeta(number_mcmc, a_prior, b_prior) #flat prior

  ### Test of model vs real
  p_test <- mean(theta_post_flat<theta0)   # larger is higher failure

  ### Number of effective sample size given shape and scale loss function
  if (two_side == 0) {
    alpha_loss <- pweibull(p_test, shape = weibull_shape, scale = weibull_scale)*alpha_max
  } else if (two_side == 1){
    p_test1    <- ifelse(p_test > 0.5, 1 - p_test, p_test)
    alpha_loss <- pweibull(p_test1, shape = weibull_shape, scale = weibull_scale)*alpha_max
  }

  return(list(alpha_loss      = alpha_loss,
              pvalue          = p_test,
              theta_post_flat = theta_post_flat,
              theta0          = theta0))
}


################################################################################
# Calculates posterior estimation for Binomial distribution given alpha_loss   #
################################################################################
theta_post_aug_bin <- function(y, N, y0, N0, alpha_loss, a0, b0,
                               number_mcmc){

  effective_N0 <- N0*alpha_loss

  a_prior  <- (y0/N0)*effective_N0 + a0
  b_prior  <- effective_N0 - (y0/N0)*effective_N0 + b0

  a_post_aug  <- y + a_prior
  b_post_aug  <- N - y + b_prior

  theta_post_aug <- rbeta(number_mcmc, a_post_aug, b_post_aug)
  return(theta_post_aug)
}


################################################################################
# Combines the loss function and posterior estimation into one function        #
################################################################################
Binomial_posterior <- function(y, N, y0, N0, alpha_max, a0, b0, number_mcmc,
                               weibull_shape, weibull_scale, two_side){

  alpha_loss         <- Loss_function(y             = y,
                                      N             = N,
                                      y0            = y0,
                                      N0            = N0,
                                      alpha_max     = alpha_max,
                                      a0            = a0,
                                      b0            = b0,
                                      number_mcmc   = number_mcmc,
                                      weibull_shape = weibull_shape,
                                      weibull_scale = weibull_scale,
                                      two_side      = two_side)

  Binomial_posterior <- theta_post_aug_bin(y           = y,
                                           N           = N,
                                           y0          = y0,
                                           N0          = N0,
                                           alpha_loss  = alpha_loss$alpha_loss,
                                           a0          = a0,
                                           b0          = b0,
                                           number_mcmc = number_mcmc)

  return(list(alpha_loss               = alpha_loss$alpha_loss,
              pvalue                   = alpha_loss$pvalue,
              Binomial_posterior       = Binomial_posterior,
              Binomial_posterior_flat  = alpha_loss$theta_post_flat,
              Binomial_prior           = alpha_loss$theta0,
              weibull_scale            = weibull_scale,
              weibull_shape            = weibull_shape,
              y                        = y,
              N                        = N,
              y0                       = y0,
              N0                       = N0,
              N0_effective             = alpha_loss$alpha_loss*N0))
}


final <- function(posterior_test){
  den_post_test     <- density(posterior_test$Binomial_posterior,adjust = 0.5)
  den_flat_test     <- density(posterior_test$Binomial_posterior_flat,adjust = 0.5)
  den_prior_test    <- density(posterior_test$Binomial_prior,adjust = 0.5)

  Testpost <- posterior_test$Binomial_posterior

  return(list(den_post_test     = den_post_test,
              den_flat_test     = den_flat_test,
              den_prior_test    = den_prior_test,
              Testpost          =  Testpost))
}




results <- function(f, posterior_test, H0, two_side, inequality){

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

  if (inequality == "<") {
    hypothesis <- paste("\"We can define mu as the mean for the test", "\n",
                        "Null Hypothesis (H_0): mu>", H0, "\n",
                        "Alternative Hypothesis (H_a): mu<", H0, "\n", "\n",
                        "P(mu<", H0, "|data)=", mean(f$Testpost < H0), "\n",
                        "We can accept H_a with a Probability of",
                        mean(f$Testpost < H0))
  }

  if (inequality == ">") {
    hypothesis <- paste("\"Define mu as the mean of the data", "\n",
                        "Null Hypothesis (H_0): mu<", H0, "\n",
                        "Alternative Hypothesis (H_a): mu>", H0, "\n", "\n",
                        "P(mu>", H0, "|data)=", mean(f$Testpost > H0), "\n",
                        "We can accept H_a with a Probability of",
                        mean(f$Testpost > H0))
  }

  ### Print
  prior_for_test_group <- list(`Effective sample size of prior (for test group)` = posterior_test$N0_effective,
                               `Bayesian p-value (new vs historical data)`       = posterior_test$pvalue,
                               `Loss function value`                             = posterior_test$alpha_loss,
                               `Sample size of prior (for test group)`           = posterior_test$N0)

  return(list(prior_for_test_group = prior_for_test_group,
              post_typeplot        = post_typeplot,
              densityplot          = densityplot,
              lossfun_plot         = lossfun_plot,
              hypothesis           = hypothesis)
  )

}




################################################################################
# Results                                                                      #
################################################################################
est <- Binomial_posterior(y, N, y0, N0, alpha_max, a0, b0, number_mcmc,
                          weibull_scale, weibull_shape, two_side)

f1 <- final(posterior_test = est)

res1 <- results(f              = f1,
                posterior_test = est,
                H0             = H0,
                two_side       = two_side,
                inequality     = inequality)


### Plot outputs
post_typeplot1 <- res1$post_typeplot
densityplot1   <- res1$densityplot
lossfun_plot1  <- res1$lossfun_plot
lossfun_plot2  <- res1$lossfun_plot


### Text outputs
hypothesis1           <- res1$hypothesis
prior_for_test_group1 <- res1$prior_for_test_group

#setClass("bdpbinomial1arm",
#         representation(post_typeplot1 = "ANY",
#                        densityplot1 = "ANY",
#                        lossfun_plot1 = "ANY",
#                        lossfun_plot2 = "ANY",
#                        hypothesis1 = "character",
#                        prior_for_test_group1 = "list"))

#me = new("bdpbinomial1arm",
#         post_typeplot1 = post_typeplot1,
#         densityplot1 = densityplot1,
#         lossfun_plot1 = lossfun_plot1,
#         lossfun_plot2 = lossfun_plot2,
#         hypothesis1 = hypothesis1,
#         prior_for_test_group1 = prior_for_test_group1)

me <- list(post_typeplot1 = post_typeplot1,
           densityplot1 = densityplot1,
           lossfun_plot1 = lossfun_plot1,
           lossfun_plot2 = lossfun_plot2,
           hypothesis1 = hypothesis1,
           prior_for_test_group1 = prior_for_test_group1)

class(me) <- "bdpbinomial1arm"

return(me)

})


#' plot
#'
#' plot
#'
#' @title plot: plot
#' @param x bdpbinomial1arm
#'
#' @examples
#'
#' @rdname plot
#' @export plot
UseMethod("plot", signature(x = "bdpbinomial1arm"), function(x){
  op <- par(ask=TRUE)
  plot(x$post_typeplot1)
  plot(x$densityplot1)
  plot(x$lossfun_plot1)
  plot(x$lossfun_plot2)
  par(op)
})

#' print
#'
#' print
#'
#' @title print: print
#' @param x bdpbinomial1arm
#'
#' @examples
#'
#' @rdname print
#' @export print
UseMethod("print", signature(x = "bdpbinomial1arm"), function(x){
  print(cat(x$hypothesis1))
  print(x$prior_for_test_group1)
})
