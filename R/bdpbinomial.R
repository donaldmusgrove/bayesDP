#' @title Bayesian Discount Prior: Binomial counts
#' @description \code{bdpbinomial} is used for estimating posterior samples from a
#'   Binomial outcome where an informative prior is used. The prior weight
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
#' @param a0 scalar. Prior value for the beta rate. Default is 1.
#' @param b0 scalar. Prior value for the beta rate. Default is 1.
#' @param two_side scalar. Indicator of two-sided test for the discount
#'   function. Default value is 1.
#'
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
#'  current control (\code{y_c} and \code{N_c}) must be input,
#'  as well as at least one of the historical treatment and historical control
#'  (\code{y0_c} and \code{N0_c}). The results
#'  are then based on the posterior distribution of the difference between
#'  current treatment and control, augmented by available historical data.
#'
#' @return \code{bdpbinomial} returns an object of class "bdpbinomial".
#' The functions \code{summary} and \code{print} are used to obtain and
#' print a summary of the results, including user inputs. The \code{plot}
#' function displays visual outputs as well.
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
#'      \item{\code{pvalue}}{
#'        numeric. The posterior probability of the stochastic comparison
#'        between the current and historical data.}
#'      \item{\code{posterior}}{
#'        vector. The posterior of the treatment group, incorporating the
#'        weighted historical data.}
#'      \item{\code{posterior_flat}}{
#'        vector. The distribution of the current treatment group, i.e., no
#'        incorporation of the historical data.}
#'      \item{\code{prior}}{
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
#' @examples
#' # One-arm trial (OPC) example
#' fit <- bdpbinomial(y_t           = 10,
#'                    N_t           = 500,
#'                    y0_t          = 25,
#'                    N0_t          = 250)
#' summary(fit)
#' print(fit)
#' #plot(fit)
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
#' #plot(fit2)
#'
#' @rdname bdpbinomial
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @aliases bdpbinomial,ANY-method
#' @export bdpbinomial
bdpbinomial <- setClass("bdpbinomial",
                        slots = c(posterior_test = "list",
                                  posterior_control = "list",
                                  f1 = "list",
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
                    a0            = 1,
                    b0            = 1,
                    number_mcmc   = 10000,
                    weibull_scale = 0.135,
                    weibull_shape = 3,
                    two_side      = 1){      #Difference margin
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
                   a0            = 1,
                   b0            = 1,
                   number_mcmc   = 10000,
                   weibull_scale = 0.135,
                   weibull_shape = 3,
                   two_side      = 1){      #Difference margin


  ################################################################################
  # Check Input                                                                  #
  ################################################################################

  intent <- c()
  if(length(y_t + N_t) != 0){
    intent <- c(intent,"current treatment")
    cat("Current Treatment\n")
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
    cat("Historical Treatment\n")
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
    cat("Current Control\n")
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
    cat("Historical Control\n")
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


  if(length(y_c + N_c + y0_c  + N0_c)!=0){
    arm2 <- TRUE
    #print("Assuming 2 arm binomial.")
  }else{
    arm2 <- FALSE
    #print("Assuming 1 arm binomial.")
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
  # Estimate weight for prior data assuming a binomial outcome                   #
  ################################################################################
  discount_function_binomial <- function(y, N, y0, N0, alpha_max, a0, b0, number_mcmc,
                                         weibull_shape, weibull_scale, two_side){

    ### Theta for using flat prior
    a_post_flat     <- y + a0
    b_post_flat     <- N - y + b0
    posterior_flat  <- rbeta(number_mcmc, a_post_flat, b_post_flat)

    ### Prior model
    a_prior  <- y0 + a0
    b_prior  <- N0 - y0 + b0
    prior    <- rbeta(number_mcmc, a_prior, b_prior) #flat prior

    ### Test of model vs real
    p_test <- mean(posterior_flat < prior)   # larger is higher failure

    ### Number of effective sample size given shape and scale discount function
    if(weibull_shape %in% c(0,Inf)){
      if(weibull_shape == 0){
        alpha_discount <- 0
      } else{
        alpha_discount <- 1
      }
    } else{
      if (two_side == 0) {
        alpha_discount <- pweibull(p_test, shape=weibull_shape, scale=weibull_scale)*alpha_max
      } else if (two_side == 1){
        p_test1    <- ifelse(p_test > 0.5, 1 - p_test, p_test)
        alpha_discount <- pweibull(p_test1, shape=weibull_shape, scale=weibull_scale)*alpha_max
      }
    }

    return(list(alpha_discount  = alpha_discount,
                pvalue          = p_test,
                posterior_flat  = posterior_flat,
                prior           = prior))
  }


  ################################################################################
  # Posterior augmentation for Binomial distribution
  ################################################################################
  posterior_augment_binomial <- function(y, N, y0, N0, alpha_discount, a0, b0,
                                         number_mcmc){

    effective_N0 <- N0 * alpha_discount

    if(is.null(N0) == FALSE){
      a_prior <- a0
      b_prior <- b0
    }else{
      a_prior <- (y0/N0)*effective_N0 + a0
      b_prior <- effective_N0 - (y0/N0)*effective_N0 + b0
    }

    a_post_aug <- y + a_prior
    b_post_aug <- N - y + b_prior

    post_aug <- rbeta(number_mcmc, a_post_aug, b_post_aug)
    return(post_aug)
  }

  ################################################################################
  # Combine discount function and posterior estimation into one function
  ################################################################################
  binomial_posterior <- function(y, N, y0, N0, alpha_max, a0, b0, number_mcmc,
                                 weibull_shape, weibull_scale, two_side){

    alpha_discount <- discount_function_binomial(y             = y,
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

    posterior <- posterior_augment_binomial(y              = y,
                                            N              = N,
                                            y0             = y0,
                                            N0             = N0,
                                            alpha_discount = alpha_discount$alpha_discount,
                                            a0             = a0,
                                            b0             = b0,
                                            number_mcmc    = number_mcmc)

    return(list(alpha_discount  = alpha_discount$alpha_discount,
                pvalue          = alpha_discount$pvalue,
                posterior       = posterior,
                posterior_flat  = alpha_discount$posterior_flat,
                prior           = alpha_discount$prior,
                weibull_scale   = weibull_scale,
                weibull_shape   = weibull_shape,
                y               = y,
                N               = N,
                y0              = y0,
                N0              = N0,
                N0_effective    = alpha_discount$alpha_discount*N0))
  }


  ################################################################################
  # Create final result class
  # - If no control, only returns posterior info for the treatment data
  ################################################################################
  final_binomial <- function(posterior_treatment, posterior_control=NULL){

    density_post_treatment  <- density(posterior_treatment$posterior,
                                       adjust = .5)
    density_flat_treatment  <- density(posterior_treatment$posterior_flat,
                                       adjust = .5)
    density_prior_treatment <- density(posterior_treatment$prior,
                                       adjust = .5)

    if(is.null(posterior_control)){
      treatment_posterior <- posterior_treatment$posterior

      return(list(density_post_treatment  = density_post_treatment,
                  density_flat_treatment  = density_flat_treatment,
                  density_prior_treatment = density_prior_treatment,
                  treatment_posterior     = treatment_posterior))
    } else{
      density_post_control  <- density(posterior_control$posterior,
                                       adjust = .5)
      density_flat_control  <- density(posterior_control$posterior_flat,
                                       adjust = .5)
      density_prior_control <- density(posterior_control$prior,
                                       adjust = .5)

      comparison_posterior <- posterior_treatment$posterior - posterior_control$posterior

      return(list(density_post_control    = density_post_control,
                  density_flat_control    = density_flat_control,
                  density_prior_control   = density_prior_control,
                  density_post_treatment  = density_post_treatment,
                  density_flat_treatment  = density_flat_treatment,
                  density_prior_treatment = density_prior_treatment,
                  comparison_posterior    = comparison_posterior))
    }

  }


  ##############################################################################
  # Run model and collect results
  ##############################################################################
  posterior_treatment <- binomial_posterior(
    y             = y_t,
    N             = N_t,
    y0            = y0_t,
    N0            = N0_t,
    alpha_max     = alpha_max[1],
    a0            = a0,
    b0            = b0,
    number_mcmc   = number_mcmc,
    weibull_scale = weibull_scale[1],
    weibull_shape = weibull_shape[1],
    two_side      = two_side)

  if (arm2==TRUE){
    posterior_control <- binomial_posterior(
      y             = y_c,
      N             = N_c,
      y0            = y0_c,
      N0            = N0_c,
      alpha_max     = alpha_max[2],
      a0            = a0,
      b0            = b0,
      number_mcmc   = number_mcmc,
      weibull_scale = weibull_scale[2],
      weibull_shape = weibull_shape[2],
      two_side      = two_side)
  }

  if (arm2==TRUE){
    f1 <- final_binomial(posterior_treatment = posterior_treatment,
                         posterior_control   = posterior_control)
  }
  else{
    f1 <- final_binomial(posterior_treatment = posterior_treatment,
                         posterior_control   = NULL)
  }

  args1 <- list(y_t           = y_t,
                N_t           = N_t,
                y0_t          = y0_t,
                N0_t          = N0_t,
                y_c           = y_c,
                N_c           = N_c,
                y0_c          = y0_c,
                N0_c          = N0_c,
                alpha_max     = alpha_max[1],
                a0            = a0,
                b0            = b0,
                number_mcmc   = number_mcmc,
                weibull_scale = weibull_scale,
                weibull_shape = weibull_shape,
                two_side      = two_side,
                arm2          = arm2,
                intent        = paste(intent,collapse=", "))

  if(arm2==TRUE){
    me <- list(posterior_treatment = posterior_treatment,
               posterior_control   = posterior_control,
               f1                  = f1,
               args1               = args1)
  } else{
      me <- list(posterior_treatment = posterior_treatment,
                 f1                  = f1,
                 args1               = args1)
  }

  class(me) <- "bdpbinomial"

  return(me)

})
