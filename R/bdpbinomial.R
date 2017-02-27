#' bdpbinomial
#' bdpbinomial
#' @title bdpbinomial: bdpbinomial
#' @param y_t Number of events for the current treatment group.
#' @param N_t Sample size of the current treatment group.
#' @param y0_t Number of events for the historical treatment group.
#' @param N0_t Sample size of the historical treatment group.
#' @param y_c Number of events for the current control group.
#' @param N_c Sample size of the current control group.
#' @param y0_c Number of events for the historical control group.
#' @param N0_c Sample size of the historical control group.
#' @param alpha_max Maximum weight the discount function can apply. Default is 1. For type="2arm", users may specify a vector of two values where the first value is used to weight the historical treatment group and the second value is used to weight the historical control group.
#' @param a0 Prior value for the beta rate. Default is 1.
#' @param b0 Prior value for the beta rate. Default is 1.
#' @param number_mcmc Number of Markov Chain Monte Carlo (MCMC) simulations. Default is 1e4.
#' @param weibull_shape Shape parameter of the Weibull discount function used to compute alpha, the weight parameter of the historical data. Default value is 3. For type="2arm", users may specify a vector of two values where the first value is used to estimate the weight of the historical treatment group and the second value is used to estimate the weight of the historical control group.
#' @param weibull_scale Scale parameter of the Weibull discount function used to compute alpha, the weight parameter of the historical data. Default value is 0.135. Two values have special treatment: 0 and Inf. For weibull_scale = 0, alpha is set to 0, i.e., no weight. For weibull_scale = Inf, alpha is set to 1, i.e., full weight. For type="2arm", users may specify a vector of two values where the first value is used to estimate the weight of the historical treatment group and the second value is used to estimate the weight of the historical control group.
#' @param two_side Indicator of two-sided test for the discount function. Default value is 1.
#' @examples
#' ### OPC (1arm) example
#' fit <- bdpbinomial(y_t           = 10,
#'                    N_t           = 500,
#'                    y0_t          = 25,
#'                    N0_t          = 250,
#'                    alpha_max     = 1,
#'                    a0            = 1,
#'                    b0            = 1,
#'                    number_mcmc   = 10000,
#'                    weibull_scale = 0.135,
#'                    weibull_shape = 3,
#'                    two_side      = 1)
#' ### Results
#' summary(fit)
#' print(fit)
#' plot(fit)
#' @rdname bdpbinomial
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
    print("Current Treatment")
  }else{
    if(is.null(y_t) == TRUE){
      print("y_t missing")
    }
    if(is.null(N_t) == TRUE){
      print("N_t missing")
    }
    stop("Current treatment not provided/incomplete.")
  }

  if(length(y0_t + N0_t) != 0){
    intent <- c(intent,"historical treatment")
    print("Historical Treatment")
  }else{
    if(length(c(y0_t, N0_t)) > 0){
      if(is.null(y0_t) == TRUE){
        print("y0_t missing")
      }
      if(is.null(N0_t) == TRUE){
        print("N0_t missing")
      }
      stop("Historical treatment incomplete.")
    }
  }

  if(length(y_c + N_c) != 0){
    intent <- c(intent,"current control")
    print("Current Control")
  }else{
    if(length(c(y_c, N_c)) > 0){
      if(is.null(y_c) == TRUE){
        print("y_c missing")
      }
      if(is.null(N_c) == TRUE){
        print("N_c missing")
      }
      stop("Current control not provided/incomplete.")
    }
  }

  if(length(y0_c + N0_c) != 0){
    intent <- c(intent,"historical control")
    print("Historical Control")
  }else{
    if(length(c(y0_c, N0_c)) > 0){
      if(is.null(y0_c) == TRUE){
        print("y0_c missing")
      }
      if(is.null(N0_c) == TRUE){
        print("N0_c missing")
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

  if (arm2==TRUE){
    f1 <- final_binomial(posterior_treatment = posterior_treatment,
                posterior_control = posterior_control)
  }
  else{
    f1 <- final_binomial(posterior_treatment = posterior_treatment,
                posterior_control = NULL)
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
                weibull_scale = weibull_scale[1],
                weibull_shape = weibull_shape[1],
                two_side      = two_side,
                arm2          = arm2,
                intent        = paste(intent,collapse=", "))

me <- list(posterior_treatment = posterior_treatment,
           posterior_control   = posterior_control,
           f1                  = f1,
           args1               = args1)

class(me) <- "bdpbinomial"

return(me)

})
