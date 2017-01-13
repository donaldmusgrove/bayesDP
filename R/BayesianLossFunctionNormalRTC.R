#' BayesianLossFunctionNormalRTC
#'
#' BayesianLossFunctionNormalRTC
#'
#' @title BayesianLossFunctionNormalRTC: BayesianLossFunctionNormalRTC
#' @param mu_t numeric
#' @param sigma2_t numeric
#' @param N_t numeric
#' @param mu_c numeric
#' @param sigma2_c numeric
#' @param N_c numeric
#' @param mu0_t numeric
#' @param sigma02_t numeric
#' @param N0_t numeric
#' @param mu0_c numeric
#' @param sigma02_c numeric
#' @param N0_c numeric
#' @param alpha_max numeric
#' @param weibull_scale numeric
#' @param weibull_shape numeric
#' @param number_mcmc numeric
#' @param two_side character
#' @param inequality character
#' @param delta numeric
#'
#' @examples
#'
#' @rdname BayesianLossFunctionNormalRTC
# @export BayesianLossFunctionNormalRTC



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


setGeneric("BayesianLossFunctionNormalRTC",
           function(mu_t = 10,
                    sigma2_t = 2,
                    N_t = 20,
                    mu_c = 15,
                    sigma2_c = 2,
                    N_c = 20,
                    mu0_t = 12,
                    sigma02_t = 1,
                    N0_t = 10,
                    mu0_c = 12,
                    sigma02_c = 1,
                    N0_c = 10,
                    alpha_max = 1,
                    weibull_scale = 0.2,
                    weibull_shape = 2,
                    number_mcmc  = 10000,
                    two_side = 1,
                    inequality = "<",
                    delta = 0){
             standardGeneric("BayesianLossFunctionNormalRTC")
           })

setMethod("BayesianLossFunctionNormalRTC",
          signature(mu_t = "numeric"),
          function(mu_t = 10,
                   sigma2_t = 2,
                   N_t = 20,
                   mu_c = 15,
                   sigma2_c = 2,
                   N_c = 20,
                   mu0_t = 12,
                   sigma02_t = 1,
                   N0_t = 10,
                   mu0_c = 12,
                   sigma02_c = 1,
                   N0_c = 10,
                   alpha_max = 1,
                   weibull_scale = 0.2,
                   weibull_shape = 2,
                   number_mcmc  = 10000,
                   two_side = 1,
                   inequality = "<",
                   delta = 0){


################################################################################
# Produce prior data weight (scalar between 0 and 1) assuming a mu outcome     #
################################################################################
Loss_function <- function(mu, sigma2, N, mu0, sigma02, N0, alpha_max, number_mcmc,
                          weibull_shape, weibull_scale, two_side){

  ### mu for using flat prior
  sigma2_post_flat <- rinvgamma(number_mcmc, (N - 1)/2, ((N - 1) * sigma2)/2)
  mu_post_flat     <- rnorm(number_mcmc, mu, (sigma2_post_flat/((N-1)+1))^0.5)

  ### Prior model (flat priors)
  sigma2_post_flat0 <- rinvgamma(number_mcmc, (N0-1)/2, ((N0-1)*sigma02)/2)
  mu_post_flat0     <- rnorm(number_mcmc, mu0, (sigma2_post_flat0/((N0-1)+1))^0.5)

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


################################################################################
# Estimate posterior for mu given alpha_loss value                             #
################################################################################
mu_post_aug <- function(mu, sigma2, N, mu0, sigma02, N0, alpha_loss,
                        number_mcmc) {
  if (N0 != 0){
    effective_N0 <- N0 * alpha_loss
    sigma2_post  <- rinvgamma(number_mcmc, (N-1)/2, ((N-1)*sigma2)/2)
    sigma2_post0 <- rinvgamma(number_mcmc, (N0-1)/2, ((N0-1)*sigma02)/2)

    mu1 <- (sigma2_post0*N*mu + sigma2_post*effective_N0*mu0)/(N*sigma2_post0 +
                                                                 sigma2_post*effective_N0)
    var_mu <- (sigma2_post*sigma2_post0)/(N*sigma2_post0 +
                                            sigma2_post*effective_N0)
  } else {
    var_mu <- rinvgamma(number_mcmc, (N - 1)/2, ((N - 1) * sigma2)/2)
    mu1    <- mu

  }
  mu_post <- rnorm(number_mcmc, mu1, sqrt(var_mu))
  return(mu_post)
}


################################################################################
# Combine loss function and posterior estimation into one function             #
################################################################################
mu_posterior <- function(mu, sigma2, N, mu0, sigma02, N0, alpha_max, number_mcmc,
                         weibull_shape, weibull_scale, two_side) {
  if (N0 != 0) {
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
  } else {
    mu_posterior <- mu_post_aug(mu          = mu,
                                sigma2      = sigma2,
                                N           = N,
                                mu0         = mu0,
                                sigma02     = sigma02,
                                N0          = N0,
                                alpha_loss  = 0,
                                number_mcmc = number_mcmc)

    alpha_loss <- list(alpha_loss   = 0,
                       pvalue       = 0,
                       mu0          = rnorm(100),
                       mu_post_flat = mu_posterior)
  }

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

final <- function(posterior_control, posterior_test) {
  den_post_control  <- density(posterior_control$mu_posterior, adjust = 0.5)
  den_flat_control  <- density(posterior_control$mu_posterior_flat,
                               adjust = 0.5)
  den_prior_control <- density(posterior_control$mu_prior, adjust = 0.5)

  den_post_test  <- density(posterior_test$mu_posterior, adjust = 0.5)
  den_flat_test  <- density(posterior_test$mu_posterior_flat, adjust = 0.5)
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
                                 scale=posterior_test$weibull_scale)*posterior_test$N0

  Loss_function_control <- pweibull(p_value,
                                    shape=posterior_control$weibull_shape,
                                    scale=posterior_control$weibull_scale)*posterior_control$N0

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
    prior_for_test_group <- list("Sample size of prior (for test group)"          = posterior_test$N0,
                                 "Effective sample size of prior(for test group)" = posterior_test$N0_effective,
                                 "Bayesian p-value (new vs historical data)"      = posterior_test$pvalue,
                                 "Loss function value"                            = posterior_test$alpha_loss)
  }

  if(posterior_control$N0==0){
    prior_for_control_group <- "No Prior Supplied"
  } else{
    prior_for_control_group <- list("Sample size of prior (for control group)"          = posterior_control$N0,
                                    "Effective sample size of prior(for control group)" = posterior_control$N0_effective,
                                    "Bayesian p-value (new vs historical data)"         = posterior_control$pvalue,
                                    "Loss function value"                               = posterior_control$alpha_loss)
  }

  return(list(prior_for_test_group    = prior_for_test_group,
              prior_for_control_group = prior_for_control_group,
              post_typeplot           = post_typeplot,
              densityplot             = densityplot,
              lossfun_plot            = lossfun_plot,
              hypothesis              = hypothesis))
}


################################################################################
# Results                                                                      #
################################################################################
#two_side   <- 1   # 0 == 1-sided, 1 === 2-sided
#inequality <- "<" # Inequality of alternate hypothesis
#N0_t       <- 10  #Number of historical subjects in test group
#N0_c       <- 0   #Number of historical subjects in control group
#delta      <- 0   #Non-inferiority zone value

posterior_test <- mu_posterior(
  mu      = mu_t,      #mean of current treatment
  sigma2  = sigma2_t,  #variance of current treatment
  N       = N_t,       #n subjects current treatment
  mu0     = mu0_t,     #mean of historical treatment
  sigma02 = sigma02_t, #variance of historical treatment
  N0      = N0_t,      #n subjects historical treatment
  alpha_max,           #Max loss function weight
  number_mcmc,         #Number of simulations to estimate posterior and loss function
  weibull_scale,       #Loss function parameter controlling the location of a weibull function
  weibull_shape,       #Loss function parameter controlling the location of a weibull function
  two_side)            #Two or one sided hypothesis test?


posterior_control <- mu_posterior(
  mu      = mu_c,      #mean of current treatment
  sigma2  = sigma2_c,  #variance of current treatment
  N       = N_c,       #n subjects current treatment
  mu0     = mu0_c,     #mean of historical treatment
  sigma02 = sigma02_c, #variance of historical treatment
  N0      = N0_c,      #n subjects historical treatment
  alpha_max,           #Max loss function weight
  number_mcmc,         #Number of simulations to estimate posterior and loss function
  weibull_scale,       #Loss function parameter controlling the location of a weibull function
  weibull_shape,       #Loss function parameter controlling the location of a weibull function
  two_side)            #Two or one sided hypothesis test?



f1 <- final(posterior_control = posterior_control,
            posterior_test    = posterior_test)


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


me <- list(post_typeplot1,
           densityplot1,
           lossfun_plot1,
           cat(hypothesis1),
           prior_for_control_group1)

setClass("BayesianLossFunctionNormalRTC",
         representation(post_typeplot1 = "ANY",
                        densityplot1 = "ANY",
                        lossfun_plot1 = "ANY",
                        hypothesis1 = "character",
                        prior_for_control_group1 = "list"))

me = new("BayesianLossFunctionNormalRTC",
         post_typeplot1 = post_typeplot1,
         densityplot1 = densityplot1,
         lossfun_plot1 = lossfun_plot1,
         hypothesis1 = hypothesis1,
         prior_for_control_group1 = prior_for_control_group1)

return(me)

})


#' plot
#'
#' plot
#'
#' @title plot: plot
#' @param x BayesianLossFunctionNormalRTC
#'
#' @examples
#'
#' @rdname plot
#' @export plot
setMethod("plot", signature(x = "BayesianLossFunctionNormalRTC"), function(x){
  op <- par(ask=TRUE)
  plot(x@post_typeplot1)
  plot(x@densityplot1)
  plot(x@lossfun_plot1)
  plot(x@lossfun_plot2)
  par(op)
})

#' plot
#'
#' plot
#'
#' @title print: print
#' @param x BayesianLossFunctionNormalRTC
#'
#' @examples
#'
#' @rdname print
#' @export print
setMethod("print", signature(x = "BayesianLossFunctionNormalRTC"), function(x){
  print(cat(x@hypothesis1))
  print(x@prior_for_test_group1)
})
