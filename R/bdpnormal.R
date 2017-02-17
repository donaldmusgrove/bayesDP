#' @importFrom MCMCpack rinvgamma
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
NULL
#' bdpnormal
#'
#' @title bdpnormal: bdpnormal
#' @param mu_t scalar. The mean of the current treatment group.
#' @param sigma_t scalar. The standard deviation of the current treatment
#'        group.
#' @param N_t scalar. The number of observations of the current treatment
#'        group.
#' @param mu0_t scalar. The mean of the historical treatment group. Required
#'        for \code{type="1arm"} (OPC trials).
#' @param sigma0_t scalar. The standard deviation of the historical treatment
#'        group. Required for \code{type="1arm"} (OPC trials).
#' @param N0_t scalar. The number of observations of the historical treatment
#'        group. Required for \code{type="1arm"} (OPC trials).
#' @param mu_c scalar. The mean of the current control group. Required for
#'        \code{type="2arm"} (RCT trials).
#' @param sigma_c scalar. The standard deviation of the current control group.
#'        Required for \code{type="2arm"} (RCT trials).
#' @param N_c scalar. The number of observations of the current control group.
#'        Required for \code{type="2arm"} (RCT trials).
#' @param mu0_c scalar. The mean of the historical control group.
#' @param sigma0_c scalar. The standard deviation of the historical control
#'        group.
#' @param N0_c scalar. The number of observations of the historical control
#'        group.
#' @param alpha_max scalar. Maximum weight the discount function can apply.
#'        Default is 1. For \code{type="2arm"}, users may optionally specify
#'        a vector of two values where the first value is used to weight the
#'        historical treatment group and the second value is used to weight
#'        the historical control group.
#' @param weibull_shape scalar. Shape parameter of the Weibull discount
#'        function used to compute alpha, the weight parameter of the historical
#'        data. Default value is 3. For \code{type="2arm"}, users may optionally
#'        specify a vector of two values where the first value is used to
#'        estimate the weight of the historical treatment group and the second
#'        value is used to estimate the weight of the historical control group.
#' @param weibull_scale scalar. Scale parameter of the Weibull discount
#'        function used to compute alpha, the weight parameter of the historical
#'        data. Default value is 0.135. Two values have special treatment:
#'        \code{0} and \code{Inf}. For \code{weibull_scale = 0}, alpha is set
#'        to 0, i.e., no weight. For \code{weibull_scale = Inf}, alpha is set
#'        to 1, i.e., full weight. For \code{type="2arm"}, users may optionally
#'        specify a vector of two values where the first value is used to
#'        estimate the weight of the historical treatment group and the second
#'        value is used to estimate the weight of the historical control group.
#' @param number_mcmc scalar. Number of Markov Chain Monte Carlo (MCMC)
#'        simulations. Default is 1e4.
#' @param two_side scalar. Indicator of two-sided test for the discount
#'        function. Default value is 1.
#'
#' @examples
#'
#' @description
#' The pdpnormal function is used for estimating posterior samples from a
#' Gaussian outcome where an informative prior is used. The prior weight
#' is determined using a discount function. This code is modeled after
#' the methodologies developed by the MDIC working group: "Informing
#' clinical trials using bench & simulations."
#'
#' @details
#' Many, many, many details to come. In fact, the best details. Believe
#' me, I know a thing or two about building details.
#'
#' @return
#' \code{bdpnormal} returns an object of class "bdpnormal".
#'
#' The functions \code{summary} and \code{print} are used to obtain and
#' print a summary of the results, including user inputs. The \code{plot}
#' function displays visual outputs as well.
#'
#' An object of class "\code{bdpnormal} " is a list containing at least
#' the following components:
#' @param posterior_treatment a list of outputs including flat, prior, and
#'        posterior samples of the (potentially) augmented treatment group.
#' @param posterior_control a list of outputs including flat, prior, and
#'        posterior samples of the (potentially) augmented control group.
#'        Not present if \code{type="1arm}.
#' @param f1 a list of posterior densities of the treatment and control
#'        groups, where the control group densities are present only if
#'        \code{type="1arm}.
#' @param args1 list of user inputs.
#'
#'
#' @rdname bdpnormal
#' @export bdpnormal

# Depends: testthat,shiny,shinyFiles,shinythemes,shinyBS,survival,ggplot2,MCMCpack,knitr,rmarkdown,parallel,MASS,arm

bdpnormal <- setClass("bdpnormal", slots = c(posterior_treatment = "list",
                                            posterior_control = "list",
                                            f1 = "list",
                                            args1 = "list"))

setGeneric("bdpnormal",
           function(mu_t = NULL,
                    sigma_t = NULL,
                    N_t = NULL,
                    mu0_t = NULL,
                    sigma0_t = NULL,
                    N0_t = NULL,
                    mu_c = NULL,
                    sigma_c = NULL,
                    N_c = NULL,
                    mu0_c = NULL,
                    sigma0_c = NULL,
                    N0_c = NULL,
                    alpha_max = 1,
                    weibull_scale = 0.135,
                    weibull_shape = 3,
                    number_mcmc  = 10000,
                    two_side = 1){
             standardGeneric("bdpnormal")
           })

setMethod("bdpnormal",
          signature(),
          function(mu_t = NULL,
                   sigma_t = NULL,
                   N_t = NULL,
                   mu0_t = NULL,
                   sigma0_t = NULL,
                   N0_t = NULL,
                   mu_c = NULL,
                   sigma_c = NULL,
                   N_c = NULL,
                   mu0_c = NULL,
                   sigma0_c = NULL,
                   N0_c = NULL,
                   alpha_max = 1,
                   weibull_scale = 0.135,
                   weibull_shape = 3,
                   number_mcmc  = 10000,
                   two_side = 1){


  ################################################################################
  # Check Input                                                                  #
  ################################################################################

  if(length(mu_c + sigma_c + N_c + mu0_c  + sigma0_c + N0_c)!=0){
    arm2 <- TRUE
    #print("Assuming 2 arm normal.")
  }else{
    arm2 <- FALSE
    #print("Assuming 1 arm normal.")
  }

  intent <- c()
  if(length(mu_t + sigma_t + N_t) != 0){
    intent <- c(intent,"current treatment")
    print("Current Treatment")
  }else{
    if(is.null(mu_t) == TRUE){
      print("mu_t missing")
    }
    if(is.null(sigma_t) == TRUE){
      print("sigma_t missing")
    }
    if(is.null(N_t) == TRUE){
      print("N_t missing")
    }
    stop("Current treatment not provided/incomplete.")
  }

  if(length(mu0_t + sigma0_t + N0_t) != 0){
    intent <- c(intent,"historical treatment")
    print("Historical Treatment")
  }else{
    if(length(c(mu0_t, sigma0_t, N0_t)) > 0){
      if(is.null(mu0_t) == TRUE){
        print("mu0_t missing")
      }
      if(is.null(sigma0_t) == TRUE){
        print("sigma0_t missing")
      }
      if(is.null(N0_t) == TRUE){
        print("N0_t missing")
      }
      stop("Historical treatment incomplete.")
    }
  }

  if(length(mu_c + sigma_c + N_c) != 0){
    intent <- c(intent,"current control")
    print("Current Control")
  }else{
    if(length(c(mu_c, sigma_c, N_c)) > 0){
      if(is.null(mu_c) == TRUE){
        print("mu_c missing")
      }
      if(is.null(sigma_c) == TRUE){
        print("sigma_c missing")
      }
      if(is.null(N_c) == TRUE){
        print("N_c missing")
      }
      stop("Current control not provided/incomplete.")
    }
  }

  if(length(mu0_c + sigma0_c + N0_c) != 0){
    intent <- c(intent,"historical control")
    print("Historical Control")
  }else{
    if(length(c(mu0_c, sigma0_c, N0_c)) > 0){
      if(is.null(mu0_c) == TRUE){
        print("mu0_c missing")
      }
      if(is.null(sigma0_c) == TRUE){
        print("sigma0_c missing")
      }
      if(is.null(N0_c) == TRUE){
        print("N0_c missing")
      }
      stop("Historical Control not provided/incomplete.")
    }
  }

  ##############################################################################
  # Quick check, if alpha_max, weibull_scale, or weibull_shape have length 1,
  # repeat input twice
  ##############################################################################
  if(length(alpha_max)==1){
    alpha_max <- rep(alpha_max, 2)
  }

  if(length(weibull_scale)==1){
    weibull_scale <- rep(weibull_scale, 2)
  }

  if(length(weibull_shape)==1){
    weibull_shape <- rep(weibull_shape, 2)
  }


  ################################################################################
  # Produce prior data weight (scalar between 0 and 1) assuming a mu outcome     #
  ################################################################################
  Discount_function <- function(mu, sigma, N, mu0, sigma0, N0, alpha_max, number_mcmc,
                            weibull_shape, weibull_scale, two_side){

    ### mu for using flat prior
    sigma2_post_flat <- rinvgamma(number_mcmc, (N - 1)/2, ((N - 1) * sigma^2)/2)
    mu_post_flat     <- rnorm(number_mcmc, mu, (sigma2_post_flat/((N-1)+1))^0.5)

    ### Prior model (flat priors)
    sigma2_post_flat0 <- rinvgamma(number_mcmc, (N0-1)/2, ((N0-1)*sigma0^2)/2)
    mu_post_flat0     <- rnorm(number_mcmc, mu0, (sigma2_post_flat0/((N0-1)+1))^0.5)

    ### Test of model vs real
    p_test <- mean(mu_post_flat < mu_post_flat0)  # larger is higher failure

    ### Number of effective sample size given shape and scale discount function
    if (two_side == 0) {
      alpha_discount <- pweibull(p_test, shape = weibull_shape, scale = weibull_scale) * alpha_max
    } else if (two_side == 1){
      p_test1    <- ifelse(p_test > 0.5, 1 - p_test, p_test)
      alpha_discount <- pweibull(p_test1, shape = weibull_shape, scale = weibull_scale) * alpha_max
    }

    return(list(alpha_discount   = alpha_discount,
                pvalue       = p_test,
                mu_post_flat = mu_post_flat,
                mu0          = mu_post_flat0))
  }


  ################################################################################
  # Estimate posterior for mu given alpha_discount value                             #
  ################################################################################
  mu_post_aug <- function(mu, sigma, N, mu0, sigma0, N0, alpha_discount,
                          number_mcmc) {
    if (is.null(N0) == FALSE){
      effective_N0 <- N0 * alpha_discount
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
  # Combine discount function and posterior estimation into one function             #
  ################################################################################
  mu_posterior <- function(mu, sigma, N, mu0, sigma0, N0, alpha_max, number_mcmc,
                           weibull_shape, weibull_scale, two_side) {
    if (is.null(N0) == FALSE){
      alpha_discount <- Discount_function(mu            = mu,
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
                                  alpha_discount  = alpha_discount$alpha_discount,
                                  number_mcmc = number_mcmc)
    } else {
      mu_posterior <- mu_post_aug(mu          = mu,
                                  sigma       = sigma,
                                  N           = N,
                                  mu0         = mu0,
                                  sigma0     = sigma0,
                                  N0          = N0,
                                  alpha_discount  = 0,
                                  number_mcmc = number_mcmc)

      alpha_discount <- list(alpha_discount   = 0,
                         pvalue       = 0,
                         mu0          = rnorm(100),
                         mu_post_flat = mu_posterior)
    }

    return(list(alpha_discount     = alpha_discount$alpha_discount,
                pvalue             = alpha_discount$pvalue,
                mu_posterior       = mu_posterior,
                mu_posterior_flat  = alpha_discount$mu_post_flat,
                mu_prior           = alpha_discount$mu0,
                weibull_scale      = weibull_scale,
                weibull_shape      = weibull_shape,
                mu                 = mu,
                N                  = N,
                mu0                = mu0,
                N0                 = N0,
                N0_effective       = alpha_discount$alpha_discount * N0))
  }

  final <- function(posterior_treatment, posterior_control = NULL) {
    if (is.null(posterior_control) == FALSE){
      den_post_control  <- density(posterior_control$mu_posterior,
                                   adjust = 0.5)
      den_flat_control  <- density(posterior_control$mu_posterior_flat,
                                   adjust = 0.5)
      den_prior_control <- density(posterior_control$mu_prior,
                                   adjust = 0.5)
    }

    den_post_treatment  <- density(posterior_treatment$mu_posterior,
                              adjust = 0.5)
    den_flat_treatment  <- density(posterior_treatment$mu_posterior_flat,
                              adjust = 0.5)
    den_prior_treatment <- density(posterior_treatment$mu_prior,
                              adjust = 0.5)

    TestMinusControl_post <- posterior_treatment$mu_posterior - posterior_control$mu_posterior
    if (is.null(N0_c) == FALSE){
    return(list(den_post_control      = den_post_control,
                den_flat_control      = den_flat_control,
                den_prior_control     = den_prior_control,
                den_post_treatment         = den_post_treatment,
                den_flat_treatment         = den_flat_treatment,
                den_prior_treatment        = den_prior_treatment,
                TestMinusControl_post = TestMinusControl_post))
    }
    else{
    return(list(den_post_treatment         = den_post_treatment,
                den_flat_treatment         = den_flat_treatment,
                den_prior_treatment        = den_prior_treatment,
                TestMinusControl_post = TestMinusControl_post))
    }
  }

  ################################################################################
  # Results                                                                      #
  ################################################################################

  posterior_treatment <- mu_posterior(
    mu      = mu_t,      #mean of current treatment
    sigma   = sigma_t,   #sd of current treatment
    N       = N_t,       #n subjects current treatment
    mu0     = mu0_t,     #mean of historical treatment
    sigma0  = sigma0_t,  #sd of historical treatment
    N0      = N0_t,      #n subjects historical treatment
    alpha_max,           #Max discount function weight
    number_mcmc,         #Number of simulations to estimate posterior and discount function
    weibull_scale,       #Discount function parameter controlling the location of a weibull function
    weibull_shape,       #Discount function parameter controlling the location of a weibull function
    two_side)            #Two or one sided hypothesis test?

  if (arm2 == TRUE){
    posterior_control <- mu_posterior(
      mu      = mu_c,      #mean of current treatment
      sigma   = sigma_c,   #sd of current treatment
      N       = N_c,       #n subjects current treatment
      mu0     = mu0_c,     #mean of historical treatment
      sigma0  = sigma0_c,  #sd of historical treatment
      N0      = N0_c,      #n subjects historical treatment
      alpha_max,           #Max discount function weight
      number_mcmc,         #Number of simulations to estimate posterior and discount function
      weibull_scale,       #Discount function parameter controlling the location of a weibull function
      weibull_shape,       #Discount function parameter controlling the location of a weibull function
      two_side)            #Two or one sided hypothesis test?
  }

  if (arm2 ==  TRUE){
    f1 <- final(posterior_treatment = posterior_treatment,
                posterior_control = posterior_control)
  }
  else{
    f1 <- final(posterior_treatment = posterior_treatment,
                posterior_control = NULL)
  }

  args1 <- list(mu_t          = mu_t,
                sigma_t       = sigma_t,
                N_t           = N_t,
                mu0_t         = mu0_t,
                sigma0_t      = sigma0_t,
                N0_t          = N0_t,
                mu_c          = mu_c,
                sigma_c       = sigma_c,
                N_c           = N_c,
                mu0_c         = mu0_c,
                sigma0_c      = sigma0_c,
                N0_c          = N0_c,
                alpha_max     = alpha_max[1],
                weibull_scale = weibull_scale[1],
                weibull_shape = weibull_shape[1],
                number_mcmc   = number_mcmc,
                two_side      = two_side,
                arm2          = arm2)

  if (arm2 == TRUE){
    me <- list(posterior_treatment = posterior_treatment,
               posterior_control = posterior_control,
               f1 = f1,
               args1 = args1)
  }
  else{
    me <- list(posterior_treatment = posterior_treatment,
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
  posterior_treatment <- x$posterior_treatment
  posterior_control <- x$posterior_control
  two_side <- x$args1$two_sid
  N0_t <- x$args1$N0_t
  N0_c <- x$args1$N0_c
  arm2 <- x$args1$arm2
  if (arm2 == TRUE){
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
  }

  D4 <- data.frame(information_sources='Posterior',
                   group="Test",
                   y=f$den_post_treatment$y,
                   x=f$den_post_treatment$x)
  D5 <- data.frame(information_sources="Current data",
                   group="Test",
                   y=f$den_flat_treatment$y,
                   x=f$den_flat_treatment$x)
  D6 <- data.frame(information_sources="Prior",
                   group="Test",
                   y=f$den_prior_treatment$y,
                   x=f$den_prior_treatment$x)

  if(is.null(N0_t) == TRUE & is.null(N0_c) == TRUE){
    D <- as.data.frame(rbind(D4,D5))
  }
  if(is.null(N0_t) == TRUE & is.null(N0_c) == FALSE){
    D <- as.data.frame(rbind(D4,D5,D1,D2,D3))
  }
  if(is.null(N0_t) == FALSE & is.null(N0_c) == TRUE){
    D <- as.data.frame(rbind(D4,D5,D6))
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

  Discount_function_treatment <- pweibull(p_value,
                                 shape=posterior_treatment$weibull_shape,
                                 scale=posterior_treatment$weibull_scale)*posterior_treatment$N0
  if(arm2 == TRUE){
    Discount_function_control <- pweibull(p_value,
                                      shape=posterior_control$weibull_shape,
                                      scale=posterior_control$weibull_scale)*posterior_control$N0
  }

  D1 <- data.frame(group="treatment",y=Discount_function_treatment,x=seq(0,1,,100))
  D2 <- data.frame(group=c("treatment"),pvalue=c(posterior_treatment$pvalue))
  D3 <- data.frame(group=c("treatment"),pvalue=c(posterior_treatment$N0_effective))

  if(arm2 == TRUE){
    D4 <- data.frame(group="control",y=Discount_function_control,x=seq(0,1,,100))
    D5 <- data.frame(group=c("control"),pvalue=c(posterior_control$pvalue))
    D6 <- data.frame(group=c("control"),pvalue=c(posterior_control$N0_effective))
  }


  discountfun_plot <- ggplot()
  if(N0_t!=0){
    discountfun_plot <- discountfun_plot +
      geom_line(data=D1,aes(y=y,x=x,colour=group),size=1) +
      geom_vline(data=D2, aes(xintercept =pvalue,colour=group),lty=2) +
      geom_hline(data=D3, aes(yintercept =pvalue,colour=group),lty=2)
  }
  if(arm2 == TRUE){
    discountfun_plot  <- discountfun_plot +
      geom_line(data=D4,aes(y=y,x=x,colour=group),size=1) +
      geom_vline(data=D5, aes(xintercept =pvalue,colour=group),lty=2) +
      geom_hline(data=D6, aes(yintercept =pvalue,colour=group),lty=2)
  }

  discountfun_plot <- discountfun_plot +
    facet_wrap(~group, ncol=1) +
    theme_bw() +
    ylab("Effective sample size for historical data") +
    xlab("Bayesian p-value (new vs historical data)")


  op <- par(ask=TRUE)
  plot(post_typeplot)
  plot(densityplot)
  plot(discountfun_plot)
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
  posterior_treatment <- x$posterior_treatment
  posterior_control <- x$posterior_control
  two_side <- x$args1$two_sid
  N0_t <- x$args1$N0_t
  N0_c <- x$args1$N0_c

  if(is.null(posterior_treatment$N0) == FALSE){
    if(posterior_treatment$N0==0){
      prior_for_treatment_group <- "No Prior Supplied"
    } else{
      prior_for_treatment_group <- list("Sample size of prior (for treatment group)"          = posterior_treatment$N0,
                                   "Effective sample size of prior(for treatment group)" = posterior_treatment$N0_effective,
                                   "Bayesian p-value (new vs historical data)"      = posterior_treatment$pvalue,
                                   "Discount function value"                            = posterior_treatment$alpha_discount)
    }
  }
  if(is.null(posterior_control$N0) == FALSE){
    if(posterior_control$N0==0){
      prior_for_control_group <- "No Prior Supplied"
    } else{
      prior_for_control_group <- list("Sample size of prior (for control group)"          = posterior_control$N0,
                                      "Effective sample size of prior(for control group)" = posterior_control$N0_effective,
                                      "Bayesian p-value (new vs historical data)"         = posterior_control$pvalue,
                                      "Discount function value"                               = posterior_control$alpha_discount)
    }
  }
  print(prior_for_treatment_group)
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
  posterior_treatment <- object$posterior_treatment
  posterior_control <- object$posterior_control
  two_side <- object$args1$two_sid
  N0_t <- object$args1$N0_t
  N0_c <- object$args1$N0_c

  if(is.null(posterior_treatment$N0) == FALSE){
    if(posterior_treatment$N0==0){
      prior_for_treatment_group <- "No Prior Supplied"
    } else{
      prior_for_treatment_group <- list("Sample size of prior (for treatment group)"          = posterior_treatment$N0,
                                   "Effective sample size of prior(for treatment group)" = posterior_treatment$N0_effective,
                                   "Bayesian p-value (new vs historical data)"      = posterior_treatment$pvalue,
                                   "Discount function value"                            = posterior_treatment$alpha_discount)
    }
  }

  if(is.null(posterior_control$N0) == FALSE){
    if(posterior_control$N0==0){
      prior_for_control_group <- "No Prior Supplied"
    } else{
      prior_for_control_group <- list("Sample size of prior (for control group)"          = posterior_control$N0,
                                      "Effective sample size of prior(for control group)" = posterior_control$N0_effective,
                                      "Bayesian p-value (new vs historical data)"         = posterior_control$pvalue,
                                      "Discount function value"                               = posterior_control$alpha_discount)
    }
  }

  print(prior_for_treatment_group)
  if(is.null(posterior_control$N0) == FALSE){
    print(prior_for_control_group)
  }

  argsdf <- suppressWarnings(data.frame(as.numeric(as.character(object$args1))))
  rownames(argsdf) <- names(object$args1)
  colnames(argsdf) <- "args"
  print(format(argsdf, scientific = FALSE))
})
