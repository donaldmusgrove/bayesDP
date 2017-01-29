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
#' bdpnormal
#'
#' @title bdpnormal: bdpnormal
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
#'
#' @examples
#'
#' @rdname bdpnormal
#' @export bdpnormal

setGeneric("bdpnormal",
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
             standardGeneric("bdpnormal")
           })

setMethod("bdpnormal",
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
  mu_post_flat     <- rnorm(number_mcmc, mu, (sigma2_post_flat/((N-1)+1))^0.5)

  ### Prior model (flat priors)
  sigma2_post_flat0 <- rinvgamma(number_mcmc, (N0-1)/2, ((N0-1)*sigma0^2)/2)
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
mu_post_aug <- function(mu, sigma, N, mu0, sigma0, N0, alpha_loss,
                        number_mcmc) {
  if (is.null(N0) == FALSE){
    effective_N0 <- N0 * alpha_loss
    sigma2_post  <- rinvgamma(number_mcmc, (N-1)/2, ((N-1)*sigma^2)/2)
    sigma2_post0 <- rinvgamma(number_mcmc, (N0-1)/2, ((N0-1)*sigma0^2)/2)

    mu1 <- (sigma2_post0*N*mu + sigma2_post*effective_N0*mu0)/(N*sigma2_post0 +
                                                                 sigma2_post*effective_N0)
    var_mu <- (sigma2_post*sigma2_post0)/(N*sigma2_post0 +
                                            sigma2_post*effective_N0)
  } else {
    var_mu <- rinvgamma(number_mcmc, (N - 1)/2, ((N - 1) * sigma^2)/2)
    mu1    <- mu

  }
  mu_post <- rnorm(number_mcmc, mu1, sqrt(var_mu))
  return(mu_post)
}


################################################################################
# Combine loss function and posterior estimation into one function             #
################################################################################
mu_posterior <- function(mu, sigma, N, mu0, sigma0, N0, alpha_max, number_mcmc,
                         weibull_shape, weibull_scale, two_side) {
  if (is.null(N0) == FALSE){
    alpha_loss <- Loss_function(mu            = mu,
                                sigma         = sigma,
                                N             = N,
                                mu0           = mu0,
                                sigma0        = sigma0,
                                N0            = N0,
                                alpha_max     = alpha_max,
                                number_mcmc   = number_mcmc,
                                weibull_shape = weibull_shape,
                                weibull_scale = weibull_scale,
                                two_side      = two_side)

    mu_posterior <- mu_post_aug(mu          = mu,
                                sigma       = sigma,
                                N           = N,
                                mu0         = mu0,
                                sigma0      = sigma0,
                                N0          = N0,
                                alpha_loss  = alpha_loss$alpha_loss,
                                number_mcmc = number_mcmc)
  } else {
    mu_posterior <- mu_post_aug(mu          = mu,
                                sigma       = sigma,
                                N           = N,
                                mu0         = mu0,
                                sigma0     = sigma0,
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
  if (is.null(posterior_control) == FALSE){
    den_post_control  <- density(posterior_control$mu_posterior,
                                 adjust = 0.5)
    den_flat_control  <- density(posterior_control$mu_posterior_flat,
                                 adjust = 0.5)
    den_prior_control <- density(posterior_control$mu_prior,
                                 adjust = 0.5)
  }

  den_post_test  <- density(posterior_test$mu_posterior,
                            adjust = 0.5)
  den_flat_test  <- density(posterior_test$mu_posterior_flat,
                            adjust = 0.5)
  den_prior_test <- density(posterior_test$mu_prior,
                            adjust = 0.5)

  TestMinusControl_post <- posterior_test$mu_posterior - posterior_control$mu_posterior
  if (is.null(N0_c) == FALSE){
  return(list(den_post_control      = den_post_control,
              den_flat_control      = den_flat_control,
              den_prior_control     = den_prior_control,
              den_post_test         = den_post_test,
              den_flat_test         = den_flat_test,
              den_prior_test        = den_prior_test,
              TestMinusControl_post = TestMinusControl_post))
  }
  else{
  return(list(den_post_test         = den_post_test,
              den_flat_test         = den_flat_test,
              den_prior_test        = den_prior_test,
              TestMinusControl_post = TestMinusControl_post))
  }
}

################################################################################
# Results                                                                      #
################################################################################

posterior_test <- mu_posterior(
  mu      = mu_t,      #mean of current treatment
  sigma   = sigma_t,   #sd of current treatment
  N       = N_t,       #n subjects current treatment
  mu0     = mu0_t,     #mean of historical treatment
  sigma0  = sigma0_t,  #sd of historical treatment
  N0      = N0_t,      #n subjects historical treatment
  alpha_max,           #Max loss function weight
  number_mcmc,         #Number of simulations to estimate posterior and loss function
  weibull_scale,       #Loss function parameter controlling the location of a weibull function
  weibull_shape,       #Loss function parameter controlling the location of a weibull function
  two_side)            #Two or one sided hypothesis test?

if (is.null(N0_c) == FALSE){
  posterior_control <- mu_posterior(
    mu      = mu_c,      #mean of current treatment
    sigma   = sigma_c,   #sd of current treatment
    N       = N_c,       #n subjects current treatment
    mu0     = mu0_c,     #mean of historical treatment
    sigma0  = sigma0_c,  #sd of historical treatment
    N0      = N0_c,      #n subjects historical treatment
    alpha_max,           #Max loss function weight
    number_mcmc,         #Number of simulations to estimate posterior and loss function
    weibull_scale,       #Loss function parameter controlling the location of a weibull function
    weibull_shape,       #Loss function parameter controlling the location of a weibull function
    two_side)            #Two or one sided hypothesis test?
}

if (is.null(N0_c) == FALSE){
  f1 <- final(posterior_control = posterior_control,
              posterior_test    = posterior_test)
}
else{
  f1 <- final(posterior_control = NULL,
              posterior_test    = posterior_test)
}

args1 <- list(mu_t = mu_t,
         sigma_t = sigma_t,
         N_t = N_t,
         mu_c = mu_c,
         sigma_c = sigma_c,
         N_c = N_c,
         mu0_t = mu0_t,
         sigma0_t = sigma0_t,
         N0_t = N0_t,
         mu0_c = mu0_c,
         sigma0_c = sigma0_c,
         N0_c = N0_c,
         alpha_max = alpha_max,
         weibull_scale = weibull_scale,
         weibull_shape = weibull_shape,
         number_mcmc  = number_mcmc,
         two_side = two_side)

if (is.null(N0_c) == FALSE){
  me <- list(posterior_test = posterior_test,
             posterior_control = posterior_control,
             f1 = f1,
             args1 = args1)
}
else{
  me <- list(posterior_test = posterior_test,
             f1 = f1,
             args1 = args1)
}

class(me) <- "bdpnormal"

return(me)

})


#' plot
#'
#' plot
#'
#' @title plot: plot
#' @param x bdpnormal
#'
#' @examples
#'
#' @rdname plot
#' @export plot
setMethod("plot", signature(x = "bdpnormal"), function(x){

  f <- x$f1
  posterior_test <- x$posterior_test
  posterior_control <- x$posterior_control
  two_side <- x$args1$two_sid
  #inequality <- x$args1$inequality
  N0_t <- x$args1$N0_t
  N0_c <- x$args1$N0_c
  #delta <- x$args1$delta

  D1 <- data.frame(information_sources='Posterior',
                   group="Control",
                   y=f$den_post_control$y,
                   x=f$den_post_control$x)
  if (is.null(N0_c) == FALSE){
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
  }
  D5 <- data.frame(information_sources="Current data",
                   group="Test",
                   y=f$den_flat_test$y,
                   x=f$den_flat_test$x)
  D6 <- data.frame(information_sources="Prior",
                   group="Test",
                   y=f$den_prior_test$y,
                   x=f$den_prior_test$x)

  if(is.null(N0_t) == TRUE & is.null(N0_c) == TRUE){
    D <- as.data.frame(rbind(D4,D5,D1,D2))
  }
  if(is.null(N0_t) == TRUE & is.null(N0_c) == FALSE){
    D <- as.data.frame(rbind(D4,D5,D1,D2,D3))
  }
  if(is.null(N0_t) == FALSE & is.null(N0_c) == TRUE){
    D <- as.data.frame(rbind(D4,D5,D6,D1,D2))
  }
  if(is.null(N0_t) == FALSE & is.null(N0_c) == FALSE){
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
#' @param x bdpnormal
#'
#' @examples
#'
#' @rdname print
#' @export print
setMethod("print", signature(x = "bdpnormal"), function(x){

  f <- x$f1
  posterior_test <- x$posterior_test
  posterior_control <- x$posterior_control
  two_side <- x$args1$two_sid
  #inequality <- x$args1$inequality
  N0_t <- x$args1$N0_t
  N0_c <- x$args1$N0_c
  #delta <- x$args1$delta

  #if(inequality=="<"){
  #  hypothesis <- paste('"We can define W as the difference between the means
  #                      for the test group versus control group, i.e. W=test
  #                      mean - control mean', '\n', 'Null Hypothesis (H_0): W>',
  #                      delta, '\n', "Alternative Hypothesis (H_a): W<", delta,
  #                      '\n', '\n', "P(W<",delta,"|data)=",
  #                      mean(f$TestMinusControl_post<delta),
  #                      '\n', "We can accept H_a with a Probability of",
  #                      mean(f$TestMinusControl_post<delta))
  #} else if(inequality==">"){
  #  hypothesis <- paste('"We can define W as the difference between the means
  #                      for the test group versus control group, i.e. W=test
  #                      mean - control mean', '\n', 'Null Hypothesis (H_0): W<',
  #                      delta, '\n', "Alternative Hypothesis (H_a): W>", delta,
  #                      '\n', '\n', "P(W<",delta,"|data)=",
  #                      mean(f$TestMinusControl_post>delta),
  #                      '\n', "We can accept H_a with a Probability of",
  #                      mean(f$TestMinusControl_post>delta))
  #}

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

  #print(cat(hypothesis))
  print(prior_for_test_group)
})

#' summary
#'
#' summary
#'
#' @title summary: summary
#' @param object bdpnormal
#'
#' @examples
#'
#' @rdname summary
#' @export summary
setMethod("summary", signature(object = "bdpnormal"), function(object){

  f <- object$f1
  posterior_test <- object$posterior_test
  posterior_control <- object$posterior_control
  two_side <- object$args1$two_sid
  #inequality <- object$args1$inequality
  N0_t <- object$args1$N0_t
  N0_c <- object$args1$N0_c
  #delta <- object$args1$delta

  #if(inequality=="<"){
  #  hypothesis <- paste('"We can define W as the difference between the means
  #                      for the test group versus control group, i.e. W=test
  #                      mean - control mean', '\n', 'Null Hypothesis (H_0): W>',
  #                      delta, '\n', "Alternative Hypothesis (H_a): W<", delta,
  #                      '\n', '\n', "P(W<",delta,"|data)=",
  #                      mean(f$TestMinusControl_post<delta),
  #                      '\n', "We can accept H_a with a Probability of",
  #                      mean(f$TestMinusControl_post<delta))
  #} else if(inequality==">"){
  #  hypothesis <- paste('"We can define W as the difference between the means
  #                      for the test group versus control group, i.e. W=test
  #                      mean - control mean', '\n', 'Null Hypothesis (H_0): W<',
  #                      delta, '\n', "Alternative Hypothesis (H_a): W>", delta,
  #                      '\n', '\n', "P(W<",delta,"|data)=",
  #                      mean(f$TestMinusControl_post>delta),
  #                      '\n', "We can accept H_a with a Probability of",
  #                      mean(f$TestMinusControl_post>delta))
  #}

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

  #print(cat(hypothesis))
  print(prior_for_test_group)
  argsdf <- data.frame(t(data.frame(object$args1)))
  names(argsdf) <- "args"
  print(argsdf)
})
