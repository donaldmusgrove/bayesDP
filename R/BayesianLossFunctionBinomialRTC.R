#' BayesianLossFunctionBinomialRTC
#'
#' BayesianLossFunctionBinomialRTC
#'
#' @title BayesianLossFunctionBinomialRTC: BayesianLossFunctionBinomialRTC
#' @param y numeric
#' @param N numeric
#' @param y0 numeric
#' @param N0 numeric
#' @param N0_max numeric
#' @param alpha_max numeric
#' @param N0_t numeric
#' @param N0_c numeric
#' @param a0 numeric
#' @param b0 numeric
#' @param number_mcmc numeric
#' @param weibull_shape numeric
#' @param weibull_scale numeric
#' @param two_side character
#' @param delta numeric
#'
#' @examples
#'
#' @rdname BayesianLossFunctionBinomialRTC
# @export BayesianLossFunctionBinomialRTC


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
#library(ggplot2)

setGeneric("BayesianLossFunctionBinomialRTC",
           function(y             = 1,          #Number of events observed  current data sets
                    N             = 400,        #Number of  current subjects
                    y0            = 10,         #Number of events observed  historical  data sets
                    N0            = 100,        #Number of historical subjects
                    N0_max        = 100,
                    alpha_max     = 1,          #Max loss function weight
                    N0_t          = 1000,
                    N0_c          = 100,
                    a0            = 1,          #Noninformative Initial priors
                    b0            = 1,          #Noninformative Initial priors
                    number_mcmc   = 10000,      #Number of simulations to estimate posterior and loss function
                    weibull_scale = .05,        #Loss function parameter controlling the location of a weibull function
                    weibull_shape = 2,          #Loss function parameter controlling the location of a weibull function
                    two_side      = 0){
             standardGeneric("BayesianLossFunctionBinomialRTC")
           })

setMethod("BayesianLossFunctionBinomialRTC",
          signature(y = "numeric"),
          function(y             = 1,          #Number of events observed  current data sets
                   N             = 400,        #Number of  current subjects
                   y0            = 10,         #Number of events observed  historical  data sets
                   N0            = 100,        #Number of historical subjects
                   N0_max        = 100,
                   alpha_max     = 1,          #Max loss function weight
                   N0_t          = 1000,
                   N0_c          = 100,
                   a0            = 1,          #Noninformative Initial priors
                   b0            = 1,          #Noninformative Initial priors
                   number_mcmc   = 10000,      #Number of simulations to estimate posterior and loss function
                   weibull_scale = .05,        #Loss function parameter controlling the location of a weibull function
                   weibull_shape = 2,          #Loss function parameter controlling the location of a weibull function
                   two_side      = 0){


################################################################################
# Produce prior data weight (scalar between 0 and 1) assuming a bin outcome    #
################################################################################
Loss_function <- function(y, N, y0, N0, alpha_max, a0, b0, number_mcmc,
                          weibull_shape, weibull_scale, two_side=0){

  ### Theta for using Flat prior
  a_post_flat     <- y + a0
  b_post_flat     <- N - y + b0
  theta_post_flat <- rbeta(number_mcmc, a_post_flat, b_post_flat)

  ### Prior model###########
  a_prior <- y0 + a0
  b_prior <- N0 - y0 + b0
  theta0  <- rbeta(number_mcmc, a_prior, b_prior) #flat prior

  ### Test of model vs real
  p_test <- mean(theta_post_flat < theta0)   #Larger is higher failure

  ### Number of effective sample size given shape and scale loss function
  if (two_side == 0) {
    alpha_loss <- pweibull(p_test, shape = weibull_shape, scale = weibull_scale) * alpha_max
  } else if (two_side == 1){
    p_test1    <- ifelse(p_test > 0.5, 1 - p_test, p_test)
    alpha_loss <- pweibull(p_test1, shape = weibull_shape, scale = weibull_scale) * alpha_max
  }

  return(list(alpha_loss      = alpha_loss,
              pvalue          = p_test,
              theta_post_flat = theta_post_flat,
              theta0          = theta0))
}


################################################################################
# Posterior estimate for Binomial distribution                                 #
################################################################################
theta_post_aug_bin <- function(y, N, y0, N0, alpha_loss, a0, b0,
                               number_mcmc){

  effective_N0 <- N0 * alpha_loss

  if(N0==0){
    a_prior <- a0
    b_prior <- b0
  }else{
    a_prior <- (y0/N0)*effective_N0 + a0
    b_prior <- effective_N0 - (y0/N0)*effective_N0 + b0
  }

  a_post_aug <- y + a_prior
  b_post_aug <- N - y + b_prior

  theta_post_aug <- rbeta(number_mcmc, a_post_aug, b_post_aug)
  return(theta_post_aug)
}


################################################################################
# Combine the loss function and posterior estimation into one function         #
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

  return(list(alpha_loss              = alpha_loss$alpha_loss,
              pvalue                  = alpha_loss$pvalue,
              Binomial_posterior      = Binomial_posterior,
              Binomial_posterior_flat = alpha_loss$theta_post_flat,
              Binomial_prior          = alpha_loss$theta0,
              weibull_scale           = weibull_scale,
              weibull_shape           = weibull_shape,
              y                       = y,
              N                       = N,
              y0                      = y0,
              N0                      = N0,
              N0_max                  = N0_max,
              N0_effective            = alpha_loss$alpha_loss * N0))
}


final <- function(posterior_control, posterior_test){
  den_post_control  <- density(posterior_control$Binomial_posterior,adjust = .5)
  den_flat_control  <- density(posterior_control$Binomial_posterior_flate,adjust = .5)
  den_prior_control <- density(posterior_control$Binomial_prior,adjust = .5)

  den_post_test     <- density(posterior_test$Binomial_posterior,adjust = .5)
  den_flat_test     <- density(posterior_test$Binomial_posterior_flate,adjust = .5)
  den_prior_test    <- density(posterior_test$Binomial_prior,adjust = .5)

  TestMinusControl_post <- posterior_test$Binomial_posterior - posterior_control$Binomial_posterior

  return(list(den_post_control      = den_post_control,
              den_flat_control      = den_flat_control,
              den_prior_control     = den_prior_control,
              den_post_test         = den_post_test,
              den_flat_test         = den_flat_test,
              den_prior_test        = den_prior_test,
              TestMinusControl_post = TestMinusControl_post))
}


results <- function(f, posterior_test, posterior_control, two_side,inequality,
                    N0_t, N0_c, delta=2){
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

### Function call & Input Variable name description
### - Input Variable name description
posterior_test <- Binomial_posterior(
  y,             #= 10,      #Number of events observed  current data sets
  N,             #= 400,     #Number of  current subjects
  y0,            #= 1,       #Number of events observed  historical  data sets
  N0 = N0_t,                 #Number of historical subjects
  alpha_max,     #= 1,       #Max loss function weight
  a0,            #= 1,       #Noninformative Initial priors
  b0,            #= 1,       #Noninformative Initial priors
  number_mcmc,   #= 10000,   #Number of simulations to estimate posterior and loss function
  weibull_scale, #= .2,      #Loss function parameter controlling the location of a weibull function
  weibull_shape, #= 2,       #Loss function parameter controlling the location of a weibull function
  two_side       #= 0
)

posterior_control <- Binomial_posterior(
  y,             #= 1,          #Number of events observed  current data sets
  N,             #= 400,        #Number of  current subjects
  y0,            #= 10,         #Number of events observed  historical  data sets
  N0 = N0_c,                    #Number of historical subjects
  alpha_max,     #= 200,        #Max loss function weight
  a0,            #= 1,          #Noninformative Initial priors
  b0,            #= 1,          #Noninformative Initial priors
  number_mcmc,   #= 10000,      #Number of simulations to estimate posterior and loss function
  weibull_scale, #= .05,        #Loss function parameter controlling the location of a weibull function
  weibull_shape, #= 2,          #Loss function parameter controlling the location of a weibull function
  two_side       #= 0
)

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

setClass("BayesianLossFunctionBinomialRTC",
         representation(post_typeplot1 = "ANY",
                        densityplot1 = "ANY",
                        lossfun_plot1 = "ANY",
                        hypothesis1 = "character",
                        prior_for_control_group1 = "list"))

me = new("BayesianLossFunctionBinomialRTC",
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
#' @param x BayesianLossFunctionBinomialRTC
#'
#' @examples
#'
#' @rdname plot
#' @export plot
setMethod("plot", signature(x = "BayesianLossFunctionBinomialRTC"), function(x){
  op <- par(ask=TRUE)
  plot(x@post_typeplot1)
  plot(x@densityplot1)
  plot(x@lossfun_plot1)
  par(op)
})

#' plot
#'
#' plot
#'
#' @title print: print
#' @param x BayesianLossFunctionBinomialRTC
#'
#' @examples
#'
#' @rdname print
#' @export print
setMethod("print", signature(x = "BayesianLossFunctionBinomialRTC"), function(x){
  print(cat(x@hypothesis1))
  print(x@prior_for_control_group1)
})
