#' bdpbinomial
#'
#' bdpbinomial
#'
#' @title bdpbinomial: bdpbinomial
#' @param y_t Number of events for the current treatment group.
#' @param N_t Sample size of the current treatment group.
#' @param y0_t Number of events for the historical treatment group.
#' @param N0_t Sample size of the historical treatment group.
#' @param y_c Number of events for the current control group.
#' @param N_c Sample size of the current control group.
#' @param y0_c Number of events for the historical control group.
#' @param N0_c Sample size of the historical control group.
#' @param type One of "1arm" or "2arm", denoting an OPC trial or a randomized control trial(RCT), respectively.
#' @param subtype subtype
#' @param alpha_max Maximum weight the discount function can apply. Default is 1. For type="2arm", users may specify a vector of two values where the first value is used to weight the historical treatment group and the second value is used to weight the historical control group.
#' @param a0 Prior value for the beta rate. Default is 1.
#' @param b0 Prior value for the beta rate. Default is 1.
#' @param number_mcmc Number of Markov Chain Monte Carlo (MCMC) simulations. Default is 1e4.
#' @param weibull_shape Shape parameter of the Weibull discount function used to compute alpha, the weight parameter of the historical data. Default value is 3. For type="2arm", users may specify a vector of two values where the first value is used to estimate the weight of the historical treatment group and the second value is used to estimate the weight of the historical control group.
#' @param weibull_scale Scale parameter of the Weibull discount function used to compute alpha, the weight parameter of the historical data. Default value is 0.135. Two values have special treatment: 0 and Inf. For weibull_scale = 0, alpha is set to 0, i.e., no weight. For weibull_scale = Inf, alpha is set to 1, i.e., full weight. For type="2arm", users may specify a vector of two values where the first value is used to estimate the weight of the historical treatment group and the second value is used to estimate the weight of the historical control group.
#' @param two_side Indicator of two-sided test for the discount function. Default value is 1.
#' @Description insert something here!
#' @Details insert something here!
#' @examples
#' ### OPC (1arm) example
#' fit <- bdpbinomial(y_t           = 10,
#'                    N_t           = 500,
#'                    y0_t          = 25,
#'                    N0_t          = 250,
#'                    type          = "1arm",
#'                    alpha_max     = 1,
#'                    a0            = 1,
#'                    b0            = 1,
#'                    number_mcmc   = 10000,
#'                    weibull_scale = 0.135,
#'                    weibull_shape = 3,
#'                    two_side      = 1)
#'
#' ### Results
#' summary(fit)
#' print(fit)
#' plot(fit)
#'
#' @rdname bdpbinomial
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
                    type = c("1arm","2arm"),
                    subtype = c(NULL,"2tc","2t","2c"),
                    alpha_max     = 1,
                    a0            = 1,
                    b0            = 1,
                    number_mcmc   = 10000,
                    weibull_scale = 0.135,
                    weibull_shape = 3,
                    two_side      = 1){
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
                   type = c("1arm","2arm"),
                   subtype = c(NULL,"2tc","2t","2c"),
                   alpha_max     = 1,
                   a0            = 1,
                   b0            = 1,
                   number_mcmc   = 10000,
                   weibull_scale = 0.135,
                   weibull_shape = 3,
                   two_side      = 1){       #Difference margin

  ################################################################################
  # Check Input                                                                  #
  ################################################################################

  if(length(y_c + N_c + y0_c  + N0_c)!=0){
    arm2 <- TRUE
    print("Assuming 2 arm.")
  }else{
    arm2 <- FALSE
    print("Assuming 1 arm.")
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

  if (arm2){
    f1 <- final_binomial(posterior_treatment = posterior_treatment,
                posterior_control = posterior_control)
  }
  else{
    f1 <- final_binomial(posterior_treatment = posterior_treatment,
                posterior_control = NULL)
  }

  args1 <- list(y_t           = y_t,
                N_t           = N_t,
                y_c           = y_c,
                N_c           = N_c,
                y0_t          = y0_t,
                N0_t          = N0_t,
                y0_c          = y0_c,
                N0_c          = N0_c,
                type          = type,
                subtype       = subtype,
                alpha_max     = alpha_max,
                a0            = a0,
                b0            = b0,
                number_mcmc   = number_mcmc,
                weibull_scale = weibull_scale,
                weibull_shape = weibull_shape,
                two_side      = two_side)

me <- list(posterior_treatment = posterior_treatment,
           posterior_control   = posterior_control,
           f1                  = f1,
           args1               = args1)

class(me) <- "bdpbinomial"

return(me)

})


#' plot
#'
#' plot
#'
#' @title plot: plot
#' @param x bdpbinomial
#'
#' @examples
#'
#' @rdname plot
#' @export plot
setMethod("plot", signature(x = "bdpbinomial"), function(x){

  f                   <- x$f1
  posterior_treatment <- x$posterior_treatment
  posterior_control   <- x$posterior_control
  two_side            <- x$args1$two_side
  N0_t                <- x$args1$N0_t
  N0_c                <- x$args1$N0_c
  arm2                <- x$args1$type
  if (arm2){
    D1 <- data.frame(information_sources='Posterior',
                     group="Control",
                     y=f$density_post_control$y,
                     x=f$density_post_control$x)

    D2 <- data.frame(information_sources="Current data",
                     group="Control",
                     y=f$density_flat_control$y,
                     x=f$density_flat_control$x)

    D3 <- data.frame(information_sources="Prior",
                     group="Control",
                     y=f$density_prior_control$y,
                     x=f$density_prior_control$x)
  }

  D4 <- data.frame(information_sources='Posterior',
                   group="Treatment",
                   y=f$density_post_treatment$y,
                   x=f$density_post_treatment$x)

  D5 <- data.frame(information_sources="Current data",
                   group="Treatment",
                   y=f$density_flat_treatment$y,
                   x=f$density_flat_treatment$x)

  D6 <- data.frame(information_sources="Prior",
                   group="Treatment",
                   y=f$density_prior_treatment$y,
                   x=f$density_prior_treatment$x)

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
    geom_line(size=2,aes(color=information_sources,lty=information_sources)) +
    theme_bw() +
    facet_wrap(~group, ncol=1,scale='free') +
    ylab("Density (PDF)") +
    xlab("values")

  densityplot <- ggplot(subset(D,information_sources=="Posterior"),
                        aes(x=x,y=y)) +
    geom_line(size=2,aes(color=group)) +
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

  discount_function_treatment <- pweibull(p_value,
    shape=posterior_treatment$weibull_shape,
    scale=posterior_treatment$weibull_scale)*posterior_treatment$N0

  discount_function_control <- pweibull(p_value,
    shape=posterior_control$weibull_shape,
    scale=posterior_control$weibull_scale)*posterior_control$N0

  D1 <- data.frame(group="Treatment",y=discount_function_treatment,x=seq(0,1,,100))
  D2 <- data.frame(group="Treatment",pvalue=c(posterior_treatment$pvalue))
  D3 <- data.frame(group="Treatment",pvalue=c(posterior_treatment$N0_effective))

  D4 <- data.frame(group="Control",y=discount_function_control,x=seq(0,1,,100))
  D5 <- data.frame(group="Control",pvalue=c(posterior_control$pvalue))
  D6 <- data.frame(group="Control",pvalue=c(posterior_control$N0_effective))


  discountfun_plot <- ggplot()
  if(N0_t!=0){
    discountfun_plot <- discountfun_plot +
      geom_line(data=D1,aes(y=y,x=x,color=group),size=1) +
      geom_vline(data=D2, aes(xintercept=pvalue,color=group),lty=2) +
      geom_hline(data=D3, aes(yintercept=pvalue,color=group),lty=2)
  }
  if(N0_c!=0){
    discountfun_plot  <- discountfun_plot +
      geom_line(data=D4,aes(y=y,x=x,color=group),size=1) +
      geom_vline(data=D5,aes(xintercept=pvalue,color=group),lty=2) +
      geom_hline(data=D6,aes(yintercept=pvalue,color=group),lty=2)
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
#' @param x bdpbinomial
#'
#' @examples
#'
#' @rdname print
#' @export print
setMethod("print", signature(x = "bdpbinomial"), function(x){

  f                   <- x$f1
  posterior_treatment <- x$posterior_treatment
  posterior_control   <- x$posterior_control
  two_side            <- x$args1$two_side
  N0_t                <- x$args1$N0_t
  N0_c                <- x$args1$N0_c


  if(is.null(posterior_treatment$N0) == FALSE){
    if(posterior_treatment$N0==0){
      prior_for_treatment_group <- "No Prior Supplied"
    }
    else{
        prior_for_treatment_group <- list(
          `Sample size of prior (for treatment group)`          = posterior_treatment$N0,
          `Effective sample size of prior(for treatment group)` = posterior_treatment$N0_effective,
          `Bayesian p-value (new vs historical data)`           = posterior_treatment$pvalue,
          `Discount function value`                             = posterior_treatment$alpha_discount)
    }
  }

  if(is.null(posterior_treatment$N0) == FALSE){
    if(posterior_treatment$N0==0){
      prior_for_control_group <- "No Prior Supplied"
    }
    else{
      prior_for_control_group <- list(
      `Sample size of prior (for control group)`          = posterior_control$N0,
      `Effective sample size of prior(for control group)` = posterior_control$N0_effective,
      `Bayesian p-value (new vs historical data)`         = posterior_control$pvalue,
      `Discount function value`                           = posterior_control$alpha_discount)
    }
  }

  print(prior_for_treatment_group)
  print(prior_for_control_group)
})

#' summary
#'
#' summary
#'
#' @title summary: summary
#' @param object bdpbinomial
#'
#' @examples
#'
#' @rdname summary
#' @export summary
setMethod("summary", signature(object = "bdpbinomial"), function(object){

  f                   <- object$f1
  posterior_treatment <- object$posterior_treatment
  posterior_control   <- object$posterior_control
  two_side            <- object$args1$two_side
  N0_t                <- object$args1$N0_t
  N0_c                <- object$args1$N0_c


  if(is.null(posterior_treatment$N0) == FALSE){
    prior_for_treatment_group <- "No Prior Supplied"
  } else{
    prior_for_treatment_group <- list(
      `Sample size of prior (for treatment group)`          = posterior_treatment$N0,
      `Effective sample size of prior(for treatment group)` = posterior_treatment$N0_effective,
      `Bayesian p-value (new vs historical data)`           = posterior_treatment$pvalue,
      `Discount function value`                             = posterior_treatment$alpha_discount)
  }

  if(is.null(posterior_treatment$N0) == FALSE){
    prior_for_control_group <- "No Prior Supplied"
  } else{
    prior_for_control_group <- list(
    `Sample size of prior (for control group)`          = posterior_control$N0,
    `Effective sample size of prior(for control group)` = posterior_control$N0_effective,
    `Bayesian p-value (new vs historical data)`         = posterior_control$pvalue,
    `Discount function value`                           = posterior_control$alpha_discount)
  }

  print(prior_for_treatment_group)
  print(prior_for_control_group)
  argsdf <- data.frame(t(data.frame(object$args1)))
  names(argsdf) <- "args"
  argsdf$"NA" <- NULL
  print(round(argsdf))
})
