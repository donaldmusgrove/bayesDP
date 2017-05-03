
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
#'   the methodologies developed in Haddad et al. (2017).
#' @param mu_t scalar. Mean of the current treatment group.
#' @param sigma_t scalar. Standard deviation of the current treatment group.
#' @param N_t scalar. Number of observations of the current treatment group.
#' @param mu0_t scalar. Mean of the historical treatment group.
#' @param sigma0_t scalar. Standard deviation of the historical treatment
#'  group.
#' @param N0_t scalar. Number of observations of the historical treatment
#'  group.
#' @param mu_c scalar. Mean of the current control group.
#' @param sigma_c scalar. Standard deviation of the current control group.
#' @param N_c scalar. Number of observations of the current control group.
#' @param mu0_c scalar. Mean of the historical control group.
#' @param sigma0_c scalar. Standard deviation of the historical control group.
#' @param N0_c scalar. Number of observations of the historical control group.
#' @param alpha_max scalar. Maximum weight the discount function can apply.
#'   Default is 1. For a two-arm trial, users may specify a vector of two values
#'   where the first value is used to weight the historical treatment group and
#'   the second value is used to weight the historical control group.
#' @param fix_alpha logical. Fix alpha at alpha_max? Default value is FALSE.
#' @param weibull_shape scalar. Shape parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 3. For a two-arm trial, users may specify a vector of two values
#'   where the first value is used to estimate the weight of the historical
#'   treatment group and the second value is used to estimate the weight of the
#'   historical control group.
#' @param weibull_scale scalar. Scale parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 0.135. For a two-arm trial, users may specify a vector of two values
#'   where the first value is used to estimate the weight of the historical
#'   treatment group and the second value is used to estimate the weight of the
#'   historical control group.
#' @param number_mcmc scalar. Number of Monte Carlo simulations. Default is 10000.
#' @param two_side logical. Indicator of two-sided test for the discount
#'   function. Default value is TRUE.
#' @details \code{bdpnormal} uses a two-stage approach for determining the
#'   strength of historical data in estimation of a mean outcome.  In the first
#'   stage, a Weibull distribution function is used as a
#'   \emph{discount function} that defines the maximum strength of the
#'   historical data (via \code{weibull_shape}, \code{weibull_scale}, and
#'   \code{alpha_max}) and discounts based on disagreement with the current data.
#'   Disagreement between current and historical data is determined by stochastically
#'   comparing the respective posterior distributions under noninformative priors.
#'   With Gaussian data, the comparison is the proability (\code{p}) that the current
#'   mean is less than the historical mean. The comparison metric \code{p} is then
#'   input into the Weibull discount function and the final strength of the
#'   historical data is returned (alpha).
#'
#'  In the second stage, posterior estimation is performed where the discount
#'  function parameter, \code{alpha}, is used as a fixed value for all posterior
#'  estimation procedures.
#'
#'  To carry out a single arm (OPC) analysis, data for the current treatment
#'  (\code{mu_t}, \code{sigma_t}, and \code{N_t}) and historical treatment
#'  (\code{mu0_t}, \code{sigma0_t}, and \code{N0_t}) must be input. The results
#'  are then based on the posterior distribution of the current data augmented
#'  by the historical data.
#'
#'  To carry our a two-arm (RCT) analysis, data for the current treatment and
#'  at least one of current or historical control data must be input.
#'  The results are then based on the posterior distribution of the difference
#'  between current treatment and control, augmented by available historical data.
#'
#'   For more details, see the \code{bdpnormal} vignette: \cr
#'   \code{vignette("bdpnormal-vignette", package="bayesDP")}
#'
#'
#' @return \code{bdpnormal} returns an object of class "bdpnormal". The
#'   functions \code{\link[=summary,bdpnormal-method]{summary}} and
#'   \code{\link[=print,bdpnormal-method]{print}} are used to obtain and print
#'   a summary of the results, including user inputs. The
#'   \code{\link[=plot,bdpnormal-method]{plot}} function displays visual
#'   outputs as well.
#'
#' An object of class \code{bdpnormal} is a list containing at least
#' the following components:
#' \describe{
#'  \item{\code{posterior_treatment}}{
#'    list. Entries contain values related to the treatment group:}
#'    \itemize{
#'      \item{\code{alpha_discount}}{
#'        numeric. Alpha value, the weighting parameter of the historical data.}
#'      \item{\code{p_hat}}{
#'        numeric. The posterior probability of the stochastic comparison
#'        between the current and historical data.}
#'      \item{\code{posterior_mu}}{
#'        vector. A vector of length \code{number_mcmc} containing the posterior
#'        mean of the treatment group. If historical treatment data is present,
#'        the posterior incorporates the weighted historical data.}
#'      \item{\code{posterior_sigma2}}{
#'        vector. A vector of length \code{number_mcmc} containing the posterior
#'        variance of the treatment group. If historical treatment data is present,
#'        the posterior incorporates the weighted historical data.}
#'      \item{\code{posterior_flat_mu}}{
#'        vector. A vector of length \code{number_mcmc} containing
#'        Monte Carlo samples of the mean of the current treatment group
#'        under a flat/non-informative prior, i.e., no incorporation of the
#'        historical data.}
#'      \item{\code{posterior_flat_sigma2}}{
#'        vector. A vector of length \code{number_mcmc} containing
#'        Monte Carlo samples of the standard deviation of the current treatment group
#'        under a flat/non-informative prior, i.e., no incorporation of the
#'        historical data.}
#'      \item{\code{prior_mu}}{
#'        vector. If historical treatment data is present, a vector of length
#'        \code{number_mcmc} containing Monte Carlo samples of the mean
#'        of the historical treatment group under a flat/non-informative prior.}
#'      \item{\code{prior_sigma2}}{
#'        vector. If historical treatment data is present, a vector of length
#'        \code{number_mcmc} containing Monte Carlo samples of the standard deviation
#'        of the historical treatment group under a flat/non-informative prior.}
#'   }
#'  \item{\code{posterior_control}}{
#'    list. Similar entries as \code{posterior_treament}. Only present if a
#'    control group is specified.
#'  }
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
#' @seealso \code{\link[=summary,bdpnormal-method]{summary}},
#'   \code{\link[=print,bdpnormal-method]{print}},
#'   and \code{\link[=plot,bdpnormal-method]{plot}} for details of each of the
#'   supported methods.
#'
#' @references
#' Haddad, T., Himes, A., Thompson, L., Irony, T., Nair, R. MDIC Computer
#'   Modeling and Simulation working group.(2017) Incorporation of stochastic
#'   engineering models as prior information in Bayesian medical device trials.
#'   \emph{Journal of Biopharmaceutical Statistics}, 1-15.
#'
#' @examples
#' # One-arm trial (OPC) example
#' fit <- bdpnormal(mu_t = 30, sigma_t = 10, N_t = 250,
#'                  mu0_t = 50, sigma0_t = 5, N0_t = 250)
#' summary(fit)
#' \dontrun{
#' plot(fit)
#' }
#'
#' # Two-arm (RCT) example
#' fit2 <- bdpnormal(mu_t = 30, sigma_t = 10, N_t = 250,
#'                   mu0_t = 50, sigma0_t = 5, N0_t = 250,
#'                   mu_c = 25, sigma_c = 10, N_c = 250,
#'                   mu0_c = 50, sigma0_c = 5, N0_c = 250)
#' summary(fit2)
#' \dontrun{
#' plot(fit2)
#' }
#'
#' @rdname bdpnormal
#' @import methods
#' @importFrom stats sd density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @aliases bdpnormal,ANY-method
#' @export bdpnormal
bdpnormal <- setClass("bdpnormal", slots = c(posterior_treatment = "list",
                                             posterior_control = "list",
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
                    fix_alpha     = FALSE,
                    weibull_scale = 0.135,
                    weibull_shape = 3,
                    number_mcmc   = 10000,
                    two_side      = TRUE){
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
                   fix_alpha     = FALSE,
                   weibull_scale = 0.135,
                   weibull_shape = 3,
                   number_mcmc   = 10000,
                   two_side      = TRUE){

  ################################################################################
  # Check Input                                                                  #
  ################################################################################

  intent <- c()
  if(length(mu_t + sigma_t + N_t) != 0){
    intent <- c(intent,"current treatment")
    #cat("Current Treatment\n")
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
    #cat("Historical Treatment\n")
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
    #cat("Current Control\n")
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
    #cat("Historical Contro\nl")
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

  if(!is.null(N_c) | !is.null(N0_c)){
    arm2 <- TRUE
  }else{
    arm2 <- FALSE
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
  # Results                                                                      #
  ################################################################################

  posterior_treatment <- posterior_normal(
    mu            = mu_t,
    sigma         = sigma_t,
    N             = N_t,
    mu0           = mu0_t,
    sigma0        = sigma0_t,
    N0            = N0_t,
    alpha_max     = alpha_max[1],
    fix_alpha     = fix_alpha,
    number_mcmc   = number_mcmc,
    weibull_scale = weibull_scale[1],
    weibull_shape = weibull_shape[1],
    two_side      = two_side)


  if (arm2){
    posterior_control <- posterior_normal(
      mu            = mu_c,
      sigma         = sigma_c,
      N             = N_c,
      mu0           = mu0_c,
      sigma0        = sigma0_c,
      N0            = N0_c,
      alpha_max     = alpha_max[2],
      fix_alpha     = fix_alpha,
      number_mcmc   = number_mcmc,
      weibull_scale = weibull_scale[2],
      weibull_shape = weibull_shape[2],
      two_side      = two_side)
  } else{
    posterior_control <- NULL
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
                alpha_max     = alpha_max,
                fix_alpha     = fix_alpha,
                weibull_scale = weibull_scale,
                weibull_shape = weibull_shape,
                number_mcmc   = number_mcmc,
                two_side      = two_side,
                arm2          = arm2,
                intent        = paste(intent,collapse=", "))

  me <- list(posterior_treatment = posterior_treatment,
             posterior_control   = posterior_control,
             args1               = args1)

  class(me) <- "bdpnormal"

  return(me)
})




################################################################################
# Normal posterior estimation
# 1) Estimate the discount function (if current+historical data both present)
# 2) Estimate the posterior of the augmented data
################################################################################
posterior_normal <- function(mu, sigma, N, mu0, sigma0, N0, alpha_max,
                             fix_alpha, number_mcmc, weibull_scale,
                             weibull_shape, two_side){

  # Compute posterior(s) of current (flat) and historical (prior) data
  # with non-informative priors
  # Current data:
  if(!is.null(N)){
    posterior_flat_sigma2 <- 1/rgamma(number_mcmc, (N - 1)/2, ((N - 1) * sigma^2)/2)
    posterior_flat_mu     <- rnorm(number_mcmc, mu, (posterior_flat_sigma2/((N-1)+1))^0.5)
  } else{
    posterior_flat_mu <- posterior_flat_sigma2 <- NULL
  }

  # Historical data:
  if(!is.null(N0)){
    prior_sigma2 <- 1/rgamma(number_mcmc, (N0-1)/2, ((N0-1)*sigma0^2)/2)
    prior_mu     <- rnorm(number_mcmc, mu0, (prior_sigma2/((N0-1)+1))^0.5)
  } else{
    prior_mu  <- prior_sigma2 <- NULL
  }

  ##############################################################################
  # Discount function
  ##############################################################################
  ### Compute stochastic comparison and alpha discount only if both
  ### N and N0 are present (i.e., current & historical data are present)
  if(!is.null(N) & !is.null(N0)){

    ### Test of model vs real
    p_hat <- mean(posterior_flat_mu < prior_mu)   # larger is higher failure

    ### Number of effective sample size given shape and scale discount function
    if(fix_alpha == TRUE){
      alpha_discount <- alpha_max
    } else{
      if (!two_side) {
        alpha_discount <- pweibull(p_hat, shape=weibull_shape,
                                   scale=weibull_scale)*alpha_max
      } else if (two_side){
        p_hat    <- ifelse(p_hat > 0.5, 1 - p_hat, p_hat)
        alpha_discount <- pweibull(p_hat, shape=weibull_shape,
                                   scale=weibull_scale)*alpha_max
      }
    }

  } else{
    alpha_discount <- NULL
    p_hat         <- NULL
  }


  ##############################################################################
  # Posterior augmentation
  # - If current or historical data are missing, this will not augment but
  #   will return the posterior of the non-missing data (with flat prior)
  ##############################################################################
  ### If only the historical data is present, compute posterior on historical
  if(is.null(N0) & !is.null(N)){
    posterior_sigma2 <- posterior_flat_sigma2
    posterior_mu     <- rnorm(number_mcmc, posterior_flat_mu, sqrt(posterior_sigma2))

  } else if(!is.null(N0) & is.null(N)){
    posterior_sigma2 <- prior_sigma2
    posterior_mu     <- rnorm(number_mcmc, prior_mu, sqrt(posterior_sigma2))

  } else if(!is.null(N0) & !is.null(N)){
    effective_N0 <- N0 * alpha_discount

    posterior_mu0 <- prior_sigma2*N*mu + posterior_flat_sigma2*effective_N0*mu0
    posterior_mu0 <- posterior_mu0 / (N*prior_sigma2 + posterior_flat_sigma2*effective_N0)

    posterior_sigma2 <- posterior_flat_sigma2*prior_sigma2
    posterior_sigma2 <- posterior_sigma2 / (N*prior_sigma2 + posterior_flat_sigma2*effective_N0)

    posterior_mu     <- rnorm(number_mcmc, posterior_mu0, sqrt(posterior_sigma2))
  }


  return(list(alpha_discount        = alpha_discount,
              p_hat                 = p_hat,
              posterior_mu          = posterior_mu,
              posterior_sigma2      = posterior_sigma2,
              posterior_flat_mu     = posterior_flat_mu,
              posterior_flat_sigma2 = posterior_flat_sigma2,
              prior_mu              = prior_mu,
              prior_sigma2          = prior_sigma2))
}

