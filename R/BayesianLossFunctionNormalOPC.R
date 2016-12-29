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
library(ggplot2)
library(MCMCpack)
library(survival)


################################################################################
### Functions
################################################################################

Loss_function <- function(mu, sigma2, N, mu0, sigma02, N0, number_mcmc,
                          weibull_shape, weibull_scale, two_side) {

  ### mu for using flat prior
  sigma2_post_flate <- rinvgamma(number_mcmc, (N-1)/2, ((N-1) * sigma2)/2)          #flat prior
  mu_post_flate     <- rnorm(number_mcmc, mu, (sigma2_post_flate/((N-1) + 1))^0.5)  #flat prior

  ### Prior model
  sigma2_post_flate0 <- rinvgamma(number_mcmc, (N0-1)/2, ((N0-1) * sigma02)/2)          #flat prior
  mu_post_flate0     <- rnorm(number_mcmc, mu0, (sigma2_post_flate0/((N0-1) + 1))^0.5)  #flat prior

  ### Test of model vs real
  p_test <- mean(mu_post_flate < mu_post_flate0)  # larger is higher failure

  ### Number of effective sample size given shape and scale loss function
  if (two_side == 0) {
    alpha_loss <- pweibull(p_test, shape = weibull_shape, scale = weibull_scale)
  } else if (two_side == 1) {
    p_test1    <- ifelse(p_test > 0.5, 1 - p_test, p_test)
    alpha_loss <- pweibull(p_test1, shape = weibull_shape, scale = weibull_scale)
  }

  return(list(alpha_loss    = alpha_loss,
              pvalue        = p_test,
              mu_post_flate = mu_post_flate,
              mu0           = mu_post_flate0))
}


################################################################################
### Calculate posterior estimation for mu distribution given alpha_loss value
### and maximum strength of (N0_max) prior if alpha_loss = 1
################################################################################
mu_post_aug <- function(mu, sigma2, N, mu0, sigma02, N0, N0_max, alpha_loss,
                        number_mcmc){
  effective_N0 <- N0_max * alpha_loss
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
mu_posterior <- function(mu, sigma2, N, mu0, sigma02, N0, N0_max, number_mcmc,
                         weibull_shape, weibull_scale, two_side) {

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
              N0_effective       = alpha_loss$alpha_loss * N0_max))
}


final <- function(posterior_control, posterior_test) {
  den_post_control  <- density(posterior_control$mu_posterior, adjust = 0.5)
  den_flat_control  <- density(posterior_control$mu_posterior_flate, adjust = 0.5)
  den_prior_control <- density(posterior_control$mu_prior, adjust = 0.5)

  den_post_test     <- density(posterior_test$mu_posterior, adjust = 0.5)
  den_flat_test     <- density(posterior_test$mu_posterior_flate, adjust = 0.5)
  den_prior_test    <- density(posterior_test$mu_prior, adjust = 0.5)

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


results <- function(f,posterior_test,H0,two_side,inequality){

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
                                 scale = posterior_test$weibull_scale)*posterior_test$N0_max

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
  prior_for_test_group <- list(`Effective sample size of prior(for test group)` = posterior_test$N0_effective,
                               `Bayesian p-value (new vs historical data)` = posterior_test$pvalue,
                               `loss function value` = posterior_test$alpha_loss,
                               N0_max = posterior_test$N0_max)

  return(list(prior_for_test_group = prior_for_test_group,
              post_typeplot        = post_typeplot,
              densityplot          = densityplot,
              lossfun_plot         = lossfun_plot,
              hypothesis           = hypothesis)
  )

}



################################################################################
### Results
################################################################################

two_side   <- 0    # 0 == 1-sided, 1 === 2-sided
H0         <- 10   # H0 value
inequality <- ">"  # Inequality of alternate hypothesis

est <- mu_posterior(mu            = 10,
                    sigma2        = 1,
                    N             = 10,   #Number of  current subjects
                    mu0           = 10,
                    sigma02       = 1,
                    N0            = 10,       #Number of historical subjects
                    N0_max        = 10,       #Maximum effective prior sample size
                    number_mcmc   = 1000,     #Number of posterior simulations
                    weibull_scale = 0.1,      #Loss function: Weibull location scale
                    weibull_shape = 1,        #Loss function: Weibull location shape
                    two_side      = two_side) # 0 == 1-sided, 1 === 2-sided

f1 <- final1(posterior_test = est)

res1 <- results(f              = f1,
                posterior_test = est,
                H0             = H0,
                two_side       = two_side,
                inequality     = inequality)


### Plot outputs
post_typeplot1 <- res1$post_typeplot

densityplot1 <- res1$densityplot

lossfun_plot1 <- res1$lossfun_plot

lossfun_plot2 <- res1$lossfun_plot

### Text outputs
hypothesis1 <- res1$hypothesis

prior_for_test_group1 <- res1$prior_for_test_group

### Display outputs
post_typeplot1
densityplot1
lossfun_plot1
lossfun_plot2

hypothesis1
prior_for_test_group1
