
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
#' @title Bayesian Discount Prior: Gaussian mean values
#' @description \code{bdpnormal} is used for estimating posterior samples from a
#'   Gaussian outcome where an informative prior is used. The prior weight
#'   is determined using a discount function. This code is modeled after
#'   the methodologies developed in Haddad (2017).
#' @param mu_t scalar. Mean of the current treatment group.
#' @param sigma_t scalar. Standard deviation of the current treatment group.
#' @param N_t scalar. Number of observations of the current treatment group.
#' @param mu0_t scalar. Mean of the historical treatment group. Required
#'   for single arm (OPC) trials.
#' @param sigma0_t scalar. Standard deviation of the historical treatment
#'  group. Required for single arm (OPC) trials.
#' @param N0_t scalar. Number of observations of the historical treatment
#'  group. Required for single arm (OPC) trials.
#' @param mu_c scalar. Mean of the current control group. Required for
#'  two arm (RCT) trials.
#' @param sigma_c scalar. Standard deviation of the current control group.
#'  Required for two arm (RCT) trials.
#' @param N_c scalar. Number of observations of the current control group.
#'  Required for two arm (RCT) trials.
#' @param mu0_c scalar. Mean of the historical control group.
#' @param sigma0_c scalar. Standard deviation of the historical control group.
#' @param N0_c scalar. Number of observations of the historical control group.
#' @param alpha_max scalar. Maximum weight the discount function can apply.
#'   Default is 1. For a two-arm trial, users may specify a vector of two values
#'   where the first value is used to weight the historical treatment group and
#'   the second value is used to weight the historical control group.
#' @param weibull_shape scalar. Shape parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 3. For a two-arm trial, users may specify a vector of two values
#'   where the first value is used to estimate the weight of the historical
#'   treatment group and the second value is used to estimate the weight of the
#'   historical control group.
#' @param weibull_scale scalar. Scale parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 0.135. Two values have special treatment: 0 and Inf. For
#'   weibull_scale = 0, alpha is set to 0, i.e., no weight. For
#'   weibull_scale = Inf, alpha is set to 1, i.e., full weight. For a two-arm
#'   trial, users may specify a vector of two values where the first value is
#'   used to estimate the weight of the historical treatment group and the
#'   second value is used to estimate the weight of the historical control
#'   group.
#' @param number_mcmc scalar. Number of Markov Chain Monte Carlo (MCMC)
#'   simulations. Default is 10000.
#' @param two_side scalar. Indicator of two-sided test for the discount
#'   function. Default value is 1.
#'
#' @details Many, many, many details to come. In fact, the best details. Believe
#' me, I know a thing or two about building details.
#'
#' @return \code{bdpnormal} returns an object of class "bdpnormal".
#' The functions \code{\link{summary}} and \code{\link{print}} are used to obtain and
#' print a summary of the results, including user inputs. The \code{\link{plot}}
#' function displays visual outputs as well.
#'
#' An object of class "\code{bdpnormal} " is a list containing at least
#' the following components:
#' \describe{
#'  \item{\code{posterior_treatment}}{
#'    list. Entries contain values related to the treatment group:}
#'    \itemize{
#'      \item{\code{alpha_discount}}{
#'        numeric. Alpha value, the weighting parameter of the historical data.}
#'      \item{\code{pvalue}}{
#'        numeric. The posterior probability of the stochastic comparison
#'        between the current and historical data.}
#'      \item{\code{mu_posterior}}{
#'        vector. The posterior of the treatment group, incorporating the
#'        weighted historical data.}
#'      \item{\code{mu_posterior_flat}}{
#'        vector. The distribution of the current treatment group, i.e., no
#'        incorporation of the historical data.}
#'      \item{\code{mu_prior}}{
#'        vector. The distribution of the historical treatment group.}
#'   }
#'  \item{\code{posterior_control}}{
#'    list. Similar entries as \code{posterior_treament}. Only present if
#'    control group is specified.}
#'  \item{\code{f1}}{
#'    list. Entries contain values related to the posterior effect:}
#'    \itemize{
#'      \item{\code{density_post_treatment}}{
#'        object of class \code{density}. Used internally to plot the density of
#'        the treatment group posterior.}
#'      \item{\code{density_flat_treatment}}{
#'        object of class \code{density}. Used internally to plot the density of
#'        the treatment group "flat" distribution.}
#'      \item{\code{density_prior_treatment}}{
#'        object of class \code{density}. Used internally to plot the density of
#'        the treatment group prior.}
#'      \item{\code{density_post_control}}{
#'        object of class \code{density}. Used internally to plot the density of
#'        the control group (if present) posterior.}
#'      \item{\code{density_flat_control}}{
#'        object of class \code{density}. Used internally to plot the density of
#'        the control group (if present) "flat" distribution.}
#'      \item{\code{density_prior_control}}{
#'        object of class \code{density}. Used internally to plot the density of
#'        the control group (if present) prior.}
#'      \item{\code{TestMinusControl_post}}{
#'        vector. If control group is present, vector contains posterior
#'        distribution of the effect estimate of treatment vs. control.
#'        control groups.}
#'   }
#'  \item{\code{args1}}{
#'    list. Entries contain user inputs. In addition, the following elements
#'    are ouput:}
#'    \itemize{
#'      \item{\code{arm2}}{
#'        binary indicator. Used internally to indicate one-arm or two-arm
#'        analysis.}
#'      \item{\code{intent}}{
#'        character. Denotes current/historical status of treatment and
#'        control groups.}
#'   }
#' }
#'
#' @references
#' Haddad, T. (2017) Incorporation of stochastic engineering models as prior
#'   information in Bayesian medical device trials.
#'   \emph{Journal of Biopharmaceutical Statistics}. To Appear.
#'
#' @examples
#' # One-arm trial (OPC) example
#' fit <- bdpnormal(mu_t = 30, sigma_t = 10, N_t = 250,
#'                  mu0_t = 50, sigma0_t = 5, N0_t = 250)
#' summary(fit)
#' plot(fit)
#'
#' # Two-arm (RCT) example
#' fit2 <- bdpnormal(mu_t = 30, sigma_t = 10, N_t = 250,
#'                   mu0_t = 50, sigma0_t = 5, N0_t = 250,
#'                   mu_c = 25, sigma_c = 10, N_c = 250,
#'                   mu0_c = 50, sigma0_c = 5, N0_c = 250)
#' summary(fit2)
#' plot(fit2)
#'
#' @rdname bdpnormal
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @aliases bdpnormal,ANY-method
#' @export bdpnormal
bdpnormal <- setClass("bdpnormal", slots = c(posterior_treatment = "list",
                                            posterior_control = "list",
                                            f1 = "list",
                                            args1 = "list"))

setGeneric("bdpnormal",
           function(mu_t          = NULL,
                    sigma_t       = NULL,
                    N_t           = NULL,
                    mu0_t         = NULL,
                    sigma0_t      = NULL,
                    N0_t          = NULL,
                    mu_c          = NULL,
                    sigma_c       = NULL,
                    N_c           = NULL,
                    mu0_c         = NULL,
                    sigma0_c      = NULL,
                    N0_c          = NULL,
                    alpha_max     = 1,
                    weibull_scale = 0.135,
                    weibull_shape = 3,
                    number_mcmc   = 10000,
                    two_side      = 1){
             standardGeneric("bdpnormal")
           })

setMethod("bdpnormal",
          signature(),
          function(mu_t          = NULL,
                   sigma_t       = NULL,
                   N_t           = NULL,
                   mu0_t         = NULL,
                   sigma0_t      = NULL,
                   N0_t          = NULL,
                   mu_c          = NULL,
                   sigma_c       = NULL,
                   N_c           = NULL,
                   mu0_c         = NULL,
                   sigma0_c      = NULL,
                   N0_c          = NULL,
                   alpha_max     = 1,
                   weibull_scale = 0.135,
                   weibull_shape = 3,
                   number_mcmc   = 10000,
                   two_side      = 1){


  ################################################################################
  # Check Input                                                                  #
  ################################################################################

  intent <- c()
  if(length(mu_t + sigma_t + N_t) != 0){
    intent <- c(intent,"current treatment")
    cat("Current Treatment\n")
  }else{
    if(is.null(mu_t) == TRUE){
      cat("mu_t missing\n")
    }
    if(is.null(sigma_t) == TRUE){
      cat("sigma_t missing\n")
    }
    if(is.null(N_t) == TRUE){
      cat("N_t missing\n")
    }
    stop("Current treatment not provided/incomplete.")
  }

  if(length(mu0_t + sigma0_t + N0_t) != 0){
    intent <- c(intent,"historical treatment")
    cat("Historical Treatment\n")
  }else{
    if(length(c(mu0_t, sigma0_t, N0_t)) > 0){
      if(is.null(mu0_t) == TRUE){
        cat("mu0_t missing\n")
      }
      if(is.null(sigma0_t) == TRUE){
        cat("sigma0_t missing\n")
      }
      if(is.null(N0_t) == TRUE){
        cat("N0_t missing\n")
      }
      stop("Historical treatment incomplete.")
    }
  }

  if(length(mu_c + sigma_c + N_c) != 0){
    intent <- c(intent,"current control")
    cat("Current Control\n")
  }else{
    if(length(c(mu_c, sigma_c, N_c)) > 0){
      if(is.null(mu_c) == TRUE){
        cat("mu_c missing\n")
      }
      if(is.null(sigma_c) == TRUE){
        cat("sigma_c missing\n")
      }
      if(is.null(N_c) == TRUE){
        cat("N_c missing\n")
      }
      stop("Current control not provided/incomplete.")
    }
  }

  if(length(mu0_c + sigma0_c + N0_c) != 0){
    intent <- c(intent,"historical control")
    cat("Historical Contro\nl")
  }else{
    if(length(c(mu0_c, sigma0_c, N0_c)) > 0){
      if(is.null(mu0_c) == TRUE){
        cat("mu0_c missing\n")
      }
      if(is.null(sigma0_c) == TRUE){
        cat("sigma0_c missing\n")
      }
      if(is.null(N0_c) == TRUE){
        cat("N0_c missing\n")
      }
      stop("Historical Control not provided/incomplete.")
    }
  }

  if(length(mu_c + sigma_c + N_c + mu0_c  + sigma0_c + N0_c)!=0){
    arm2 <- TRUE
    #print("Assuming 2 arm normal.")
  }else{
    arm2 <- FALSE
    #print("Assuming 1 arm normal.")
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
    sigma2_post_flat <- 1/rgamma(number_mcmc, (N - 1)/2, ((N - 1) * sigma^2)/2)
    mu_post_flat     <- rnorm(number_mcmc, mu, (sigma2_post_flat/((N-1)+1))^0.5)

    ### Prior model (flat priors)
    sigma2_post_flat0 <- 1/rgamma(number_mcmc, (N0-1)/2, ((N0-1)*sigma0^2)/2)
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
      sigma2_post  <- 1/rgamma(number_mcmc, (N-1)/2, ((N-1)*sigma^2)/2)
      sigma2_post0 <- 1/rgamma(number_mcmc, (N0-1)/2, ((N0-1)*sigma0^2)/2)

      mu1 <- (sigma2_post0*N*mu + sigma2_post*effective_N0*mu0)/(N*sigma2_post0 +
                                                                   sigma2_post*effective_N0)
      var_mu <- (sigma2_post*sigma2_post0)/(N*sigma2_post0 +
                                              sigma2_post*effective_N0)
    } else {
      var_mu <- 1/rgamma(number_mcmc, (N - 1)/2, ((N - 1) * sigma^2)/2)
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

      mu_posterior <- mu_post_aug(mu             = mu,
                                  sigma          = sigma,
                                  N              = N,
                                  mu0            = mu0,
                                  sigma0         = sigma0,
                                  N0             = N0,
                                  alpha_discount = alpha_discount$alpha_discount,
                                  number_mcmc    = number_mcmc)
    } else {
      mu_posterior <- mu_post_aug(mu             = mu,
                                  sigma          = sigma,
                                  N              = N,
                                  mu0            = mu0,
                                  sigma0         = sigma0,
                                  N0             = N0,
                                  alpha_discount = 0,
                                  number_mcmc    = number_mcmc)

      alpha_discount <- list(alpha_discount = 0,
                         pvalue             = 0,
                         mu0                = rnorm(100),
                         mu_post_flat       = mu_posterior)
    }

    return(list(alpha_discount    = alpha_discount$alpha_discount,
                pvalue            = alpha_discount$pvalue,
                mu_posterior      = mu_posterior,
                mu_posterior_flat = alpha_discount$mu_post_flat,
                mu_prior          = alpha_discount$mu0))
  }

  final <- function(posterior_treatment, posterior_control = NULL) {
    if (is.null(posterior_control) == FALSE){
      density_post_control  <- density(posterior_control$mu_posterior,
                                   adjust = 0.5)
      density_flat_control  <- density(posterior_control$mu_posterior_flat,
                                   adjust = 0.5)
      density_prior_control <- density(posterior_control$mu_prior,
                                   adjust = 0.5)
    }

    density_post_treatment  <- density(posterior_treatment$mu_posterior,
                              adjust = 0.5)
    density_flat_treatment  <- density(posterior_treatment$mu_posterior_flat,
                              adjust = 0.5)
    density_prior_treatment <- density(posterior_treatment$mu_prior,
                              adjust = 0.5)

    TestMinusControl_post <- posterior_treatment$mu_posterior - posterior_control$mu_posterior

    if (is.null(N0_c) == FALSE){
    return(list(density_post_control    = density_post_control,
                density_flat_control    = density_flat_control,
                density_prior_control   = density_prior_control,
                density_post_treatment  = density_post_treatment,
                density_flat_treatment  = density_flat_treatment,
                density_prior_treatment = density_prior_treatment,
                TestMinusControl_post   = TestMinusControl_post))
    }
    else{
    return(list(density_post_treatment    = density_post_treatment,
                density_flat_treatment    = density_flat_treatment,
                density_prior_treatment   = density_prior_treatment,
                TestMinusControl_post = TestMinusControl_post))
    }
  }

  ################################################################################
  # Results                                                                      #
  ################################################################################

  posterior_treatment <- mu_posterior(
    mu            = mu_t,
    sigma         = sigma_t,
    N             = N_t,
    mu0           = mu0_t,
    sigma0        = sigma0_t,
    N0            = N0_t,
    alpha_max     = alpha_max[1],
    number_mcmc   = number_mcmc,
    weibull_scale = weibull_scale[1],
    weibull_shape = weibull_shape[1],
    two_side      = two_side)


  if (arm2 == TRUE){
    posterior_control <- mu_posterior(
      mu            = mu_c,
      sigma         = sigma_c,
      N             = N_c,
      mu0           = mu0_c,
      sigma0        = sigma0_c,
      N0            = N0_c,
      alpha_max     = alpha_max[2],
      number_mcmc   = number_mcmc,
      weibull_scale = weibull_scale[2],
      weibull_shape = weibull_shape[2],
      two_side      = two_side)
  }

  if (arm2 ==  TRUE){
    f1 <- final(posterior_treatment = posterior_treatment,
                posterior_control = posterior_control)
  }
  else{
    f1 <- final(posterior_treatment = posterior_treatment,
                posterior_control   = NULL)
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
                weibull_scale = weibull_scale,
                weibull_shape = weibull_shape,
                number_mcmc   = number_mcmc,
                two_side      = two_side,
                arm2          = arm2,
                intent        = paste(intent,collapse=", "))

  if (arm2 == TRUE){
    me <- list(posterior_treatment = posterior_treatment,
               posterior_control   = posterior_control,
               f1                  = f1,
               args1               = args1)
  }
  else{
    me <- list(posterior_treatment = posterior_treatment,
               f1                  = f1,
               args1               = args1)
  }

  class(me) <- "bdpnormal"

  return(me)

})
