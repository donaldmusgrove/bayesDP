#' @title Bayesian Discount Prior: Binomial counts
#' @description \code{bdpbinomial} is used for estimating posterior samples from a
#'   binomial outcome where an informative prior is used. The prior weight
#'   is determined using a discount function. This code is modeled after
#'   the methodologies developed in Haddad et al. (2017).
#' @param y_t scalar. Number of events for the current treatment group.
#' @param N_t scalar. Sample size of the current treatment group.
#' @param y0_t scalar. Number of events for the historical treatment group.
#' @param N0_t scalar. Sample size of the historical treatment group.
#' @param y_c scalar. Number of events for the current control group.
#' @param N_c scalar. Sample size of the current control group.
#' @param y0_c scalar. Number of events for the historical control group.
#' @param N0_c scalar. Sample size of the historical control group.
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
#' @param a0 scalar. Prior value for the beta rate. Default is 1.
#' @param b0 scalar. Prior value for the beta rate. Default is 1.
#' @param two_side logical. Indicator of two-sided test for the discount
#'   function. Default value is TRUE.
#' @details \code{bdpbinomial} uses a two-stage approach for determining the
#'   strength of historical data in estimation of a binomial count mean outcome.
#'   In the first stage, a Weibull distribution function is used as a
#'   \emph{discount function} that defines the maximum strength of the
#'   historical data (via \code{weibull_shape}, \code{weibull_scale}, and
#'   \code{alpha_max}) and discounts based on disagreement with the current data.
#'   Disagreement between current and historical data is determined by stochastically
#'   comparing the respective posterior distributions under noninformative priors.
#'   With binomial data, the comparison is the proability (\code{p}) that the current
#'   count is less than the historical count. The comparison metric \code{p} is then
#'   input into the Weibull discount function and the final strength of the
#'   historical data is returned (alpha).
#'
#'  In the second stage, posterior estimation is performed where the discount
#'  function parameter, \code{alpha}, is used as a fixed value for all posterior
#'  estimation procedures.
#'
#'  To carry out a single arm (OPC) analysis, data for the current treatment
#'  (\code{y_t} and \code{N_t}) and historical treatment (\code{y0_t} and
#'  \code{N0_t}) must be input. The results are then based on the posterior
#'  distribution of the current data augmented by the historical data.
#'
#'  To carry our a two-arm (RCT) analysis, data for the current treatment and
#'  at least one of current or historical control data must be input. The results
#'  are then based on the posterior distribution of the difference between
#'  current treatment and control, augmented by available historical data.
#'
#'   For more details, see the \code{bdpbinomial} vignette: \cr
#'   \code{vignette("bdpbinomial-vignette", package="bayesDP")}
#'
#' @return \code{bdpbinomial} returns an object of class "bdpbinomial". The
#'   functions \code{\link[=summary,bdpbinomial-method]{summary}} and
#'   \code{\link[=print,bdpbinomial-method]{print}} are used to obtain and
#'   print a summary of the results, including user inputs. The
#'   \code{\link[=plot,bdpbinomial-method]{plot}} function displays visual
#'   outputs as well.
#'
#' An object of class \code{bdpbinomial} is a list containing at least
#' the following components:
#'
#' \describe{
#'  \item{\code{posterior_treatment}}{
#'    list. Entries contain values related to the treatment group:}
#'    \itemize{
#'      \item{\code{alpha_discount}}{
#'        numeric. Alpha value, the weighting parameter of the historical data.}
#'      \item{\code{p_hat}}{
#'        numeric. The posterior probability of the stochastic comparison
#'        between the current and historical data.}
#'      \item{\code{posterior}}{
#'        vector. A vector of length \code{number_mcmc} containing
#'        posterior Monte Carlo samples of the event rate of the treatment
#'        group. If historical treatment data is present, the posterior
#'        incorporates the weighted historical data.}
#'      \item{\code{posterior_flat}}{
#'        vector. A vector of length \code{number_mcmc} containing
#'        Monte Carlo samples of the event rate of the current treatment group
#'        under a flat/non-informative prior, i.e., no incorporation of the
#'        historical data.}
#'      \item{\code{prior}}{
#'        vector. If historical treatment data is present, a vector of length
#'        \code{number_mcmc} containing Monte Carlo samples of the event rate
#'        of the historical treatment group under a flat/non-informative prior.}
#'   }
#'  \item{\code{posterior_control}}{
#'    list. Similar entries as \code{posterior_treament}. Only present if a
#'    control group is specified.}
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
#' @seealso \code{\link[=summary,bdpbinomial-method]{summary}},
#'   \code{\link[=print,bdpbinomial-method]{print}},
#'   and \code{\link[=plot,bdpbinomial-method]{plot}} for details of each of the
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
#' fit <- bdpbinomial(y_t           = 10,
#'                    N_t           = 500,
#'                    y0_t          = 25,
#'                    N0_t          = 250)
#' summary(fit)
#' print(fit)
#' \dontrun{
#' plot(fit)
#' }
#'
#' # Two-arm (RCT) example
#' fit2 <- bdpbinomial(y_t = 10,
#'                     N_t = 500,
#'                     y0_t = 25,
#'                     N0_t = 250,
#'                     y_c = 8,
#'                     N_c = 500,
#'                     y0_c = 20,
#'                     N0_c = 250)
#' summary(fit2)
#' print(fit2)
#' \dontrun{
#' plot(fit2)
#' }
#'
#' @rdname bdpbinomial
#' @import methods
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @aliases bdpbinomial,ANY-method
#' @export bdpbinomial
bdpbinomial <- setClass("bdpbinomial",
                        slots = c(posterior_test = "list",
                                  posterior_control = "list",
                                  args1 = "list"))

setGeneric("bdpbinomial",
           function(y_t           = NULL,
                    N_t           = NULL,
                    y0_t          = NULL,
                    N0_t          = NULL,
                    y_c           = NULL,
                    N_c           = NULL,
                    y0_c          = NULL,
                    N0_c          = NULL,
                    alpha_max     = 1,
                    fix_alpha     = FALSE,
                    a0            = 1,
                    b0            = 1,
                    number_mcmc   = 10000,
                    weibull_scale = 0.135,
                    weibull_shape = 3,
                    two_side      = TRUE){
             standardGeneric("bdpbinomial")
           })

setMethod("bdpbinomial",
          signature(),
          function(y_t           = NULL,
                   N_t           = NULL,
                   y0_t          = NULL,
                   N0_t          = NULL,
                   y_c           = NULL,
                   N_c           = NULL,
                   y0_c          = NULL,
                   N0_c          = NULL,
                   alpha_max     = 1,
                   fix_alpha     = FALSE,
                   a0            = 1,
                   b0            = 1,
                   number_mcmc   = 10000,
                   weibull_scale = 0.135,
                   weibull_shape = 3,
                   two_side      = TRUE){


  ################################################################################
  # Check Input                                                                  #
  ################################################################################

  intent <- c()
  if(length(y_t + N_t) != 0){
    intent <- c(intent,"current treatment")
    #cat("Current Treatment\n")
  }else{
    if(is.null(y_t) == TRUE){
      cat("y_t missing\n")
    }
    if(is.null(N_t) == TRUE){
      cat("N_t missing\n")
    }
    stop("Current treatment not provided/incomplete.")
  }

  if(length(y0_t + N0_t) != 0){
    intent <- c(intent,"historical treatment")
    #cat("Historical Treatment\n")
  }else{
    if(length(c(y0_t, N0_t)) > 0){
      if(is.null(y0_t) == TRUE){
        cat("y0_t missing\n")
      }
      if(is.null(N0_t) == TRUE){
        cat("N0_t missing\n")
      }
      stop("Historical treatment incomplete.")
    }
  }

  if(length(y_c + N_c) != 0){
    intent <- c(intent,"current control")
    #cat("Current Control\n")
  }else{
    if(length(c(y_c, N_c)) > 0){
      if(is.null(y_c) == TRUE){
        cat("y_c missing\n")
      }
      if(is.null(N_c) == TRUE){
        cat("N_c missing\n")
      }
      stop("Current control not provided/incomplete.")
    }
  }

  if(length(y0_c + N0_c) != 0){
    intent <- c(intent,"historical control")
    #cat("Historical Control\n")
  }else{
    if(length(c(y0_c, N0_c)) > 0){
      if(is.null(y0_c) == TRUE){
        cat("y0_c missing\n")
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




  ##############################################################################
  # Run model and collect results
  ##############################################################################
  posterior_treatment <- posterior_binomial(
    y             = y_t,
    N             = N_t,
    y0            = y0_t,
    N0            = N0_t,
    alpha_max     = alpha_max[1],
    fix_alpha     = fix_alpha,
    a0            = a0,
    b0            = b0,
    number_mcmc   = number_mcmc,
    weibull_scale = weibull_scale[1],
    weibull_shape = weibull_shape[1],
    two_side      = two_side)

  if(arm2){
    posterior_control <- posterior_binomial(
      y             = y_c,
      N             = N_c,
      y0            = y0_c,
      N0            = N0_c,
      alpha_max     = alpha_max[2],
      fix_alpha     = fix_alpha,
      a0            = a0,
      b0            = b0,
      number_mcmc   = number_mcmc,
      weibull_scale = weibull_scale[2],
      weibull_shape = weibull_shape[2],
      two_side      = two_side)
  } else{
    posterior_control <- NULL
  }


  args1 <- list(y_t           = y_t,
                N_t           = N_t,
                y0_t          = y0_t,
                N0_t          = N0_t,
                y_c           = y_c,
                N_c           = N_c,
                y0_c          = y0_c,
                N0_c          = N0_c,
                alpha_max     = alpha_max,
                fix_alpha     = fix_alpha,
                a0            = a0,
                b0            = b0,
                number_mcmc   = number_mcmc,
                weibull_scale = weibull_scale,
                weibull_shape = weibull_shape,
                two_side      = two_side,
                arm2          = arm2,
                intent        = paste(intent,collapse=", "))

  me <- list(posterior_treatment = posterior_treatment,
             posterior_control   = posterior_control,
             args1               = args1)

  class(me) <- "bdpbinomial"

  return(me)

})




################################################################################
# Binomial posterior estimation
# 1) Estimate the discount function (if current+historical data both present)
# 2) Estimate the posterior of the augmented data
################################################################################
posterior_binomial <- function(y, N, y0, N0, alpha_max, fix_alpha, a0, b0,
                               number_mcmc, weibull_scale, weibull_shape,
                               two_side){

  ### Compute posterior(s) of current (flat) and historical (prior) data
  ### with non-informative priors
  if(!is.null(N)){
    posterior_flat  <- rbeta(number_mcmc, y + a0, N - y + b0)
  } else{
    posterior_flat <- NULL
  }

  if(!is.null(N0)){
    prior    <- rbeta(number_mcmc, y0 + a0, N0 - y0 + b0)
  } else{
    prior <- NULL
  }

  ##############################################################################
  # Discount function
  ##############################################################################
  ### Compute stochastic comparison and alpha discount only if both
  ### N and N0 are present (i.e., current & historical data are present)
  if(!is.null(N) & !is.null(N0)){

    ### Test of model vs real
    p_hat <- mean(posterior_flat < prior)   # larger is higher failure

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
    posterior <- posterior_flat

  } else if(!is.null(N0) & is.null(N)){
    posterior <- prior

  } else if(!is.null(N0) & !is.null(N)){
    effective_N0 <- N0 * alpha_discount
    a_prior    <- (y0/N0)*effective_N0 + a0
    b_prior    <- effective_N0 - (y0/N0)*effective_N0 + b0
    a_post_aug <- y + a_prior
    b_post_aug <- N - y + b_prior
    posterior  <- rbeta(number_mcmc, a_post_aug, b_post_aug)
  }

  return(list(alpha_discount  = alpha_discount,
              p_hat           = p_hat,
              posterior       = posterior,
              posterior_flat  = posterior_flat,
              prior           = prior,
              weibull_scale   = weibull_scale,
              weibull_shape   = weibull_shape,
              y               = y,
              N               = N,
              y0              = y0,
              N0              = N0))
}


