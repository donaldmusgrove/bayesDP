#' This code is used for estimating posterior samples from a binary outcome
#' where an informative prior is used. The prior weight is determined using a
#' discount function. In addition this code simulate many trials in order to get
#' trial characteristics you must specify the parameters of the discount function
#' as well as the maximum strength for the prior. This code assumes a
#' non-adaptive trial.
#' This code is modeled after the methodologies developed by the MDIC working
#' group: "Informing clinical trials using bench & simulations"
#' Developer: Tarek Haddad
#' Tarek.D.Haddad@Medtronic.com
#' Last modified:1/26/2016
#'
#' bdpbinomial1arm
#'
#' @title bdpbinomial1arm: bdpbinomial1arm
#' @param y_t numeric
#' @param N_t numeric
#' @param y0_t numeric
#' @param N0_t numeric
#' @param alpha_max numeric
#' @param a0 numeric
#' @param b0 numeric
#' @param number_mcmc numeric
#' @param weibull_shape numeric
#' @param weibull_scale numeric
#' @param two_side character
#'
#' @examples
#'
#' @rdname bdpbinomial1arm
#' @export bdpbinomial1arm

setGeneric("bdpbinomial1arm",
           function(y_t           = NULL,
                    N_t           = NULL,
                    y0_t          = NULL,
                    N0_t          = NULL,
                    alpha_max     = 1,
                    a0            = 1,
                    b0            = 1,
                    number_mcmc   = 10000,
                    weibull_scale = 0.135,
                    weibull_shape = 3,
                    two_side      = 1){
             standardGeneric("bdpbinomial1arm")
           })

setMethod("bdpbinomial1arm",
          signature(),
          function(y_t           = NULL,
                   N_t           = NULL,
                   y0_t          = NULL,
                   N0_t          = NULL,
                   alpha_max     = 1,
                   a0            = 1,
                   b0            = 1,
                   number_mcmc   = 10000,
                   weibull_scale = 0.135,
                   weibull_shape = 3,
                   two_side      = 1){


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

  f1 <- final_binomial(posterior_treatment = posterior_treatment)

  args1 <- list(y_t           = y_t,
                N_t           = N_t,
                y0_t          = y0_t,
                N0_t          = N0_t,
                alpha_max     = alpha_max[1],
                a0            = a0,
                b0            = b0,
                number_mcmc   = number_mcmc,
                weibull_scale = weibull_scale[1],
                weibull_shape = weibull_shape[1],
                two_side      = two_side)

  me <- list(posterior_treatment = posterior_treatment,
             f1                  = f1,
             args1               = args1)

  class(me) <- "bdpbinomial1arm"

  return(me)

})

#' plot
#'
#' plot
#'
#' @title plot: plot
#' @param x bdpbinomial1arm
#'
#' @examples
#'
#' @method plot
#'
#' @rdname plot
#' @export plot
setMethod("plot", signature(x = "bdpbinomial1arm"), function(x){

  f                   <- x$f1
  posterior_treatment <- x$posterior_treatment
  two_side            <- x$args1$two_side

  D4 <- data.frame(information_sources = "Posterior",
                   group               = "Treatment",
                   y                   = f$density_post_treatment$y,
                   x                   = f$density_post_treatment$x)

  D5 <- data.frame(information_sources = "Current data",
                   group               = "Treatment",
                   y                   = f$density_flat_treatment$y,
                   x                   = f$density_flat_treatment$x)

  D6 <- data.frame(information_sources = "Prior",
                   group               = "Treatment",
                   y                   = f$density_prior_treatment$y,
                   x                   = f$density_prior_treatment$x)

  D <- as.data.frame(rbind(D4, D5, D6))

  D$information_sources <- factor(D$information_sources,
                                  levels = (c("Posterior", "Current data", "Prior")))

  post_typeplot <- ggplot(D, aes(x = x, y = y)) +
    geom_line(size = 2, aes(colour = information_sources, lty = information_sources)) +
    theme_bw() +
    facet_wrap(~group, ncol = 1, scale = "free") +
    ylab("Density (PDF)") +
    xlab("values")

  densityplot <- ggplot(subset(D, information_sources == "Posterior"), aes(x = x, y = y)) +
    geom_line(size = 2, aes(colour = group)) +
    ylab("Density (PDF)") +
    xlab("values") +
    theme_bw()

  if (two_side == 1) {
    p_value = seq(0, 1, , 100)
    p_value = ifelse(p_value > 0.5, 1 - p_value, p_value)
  }
  if (two_side == 0) {
    p_value = seq(0, 1, , 100)
  }

  discount_function_treatment <- pweibull(p_value,
    shape = posterior_treatment$weibull_shape,
    scale = posterior_treatment$weibull_scale)*posterior_treatment$N0

  D1 <- data.frame(group = "Treatment", y = discount_function_treatment, x = seq(0, 1, , 100))
  D2 <- data.frame(group = c("Treatment"), pvalue = c(posterior_treatment$pvalue))
  D3 <- data.frame(group = c("Treatment"), pvalue = c(posterior_treatment$N0_effective))

  discountfun_plot <- ggplot() +
    geom_line(data = D1, aes(y = y, x = x, colour = group), size = 1) +
    geom_vline(data = D2, aes(xintercept = pvalue, colour = group), lty = 2) +
    geom_hline(data = D3, aes(yintercept = pvalue, colour = group), lty = 2) +
    facet_wrap(~group, ncol = 1) +
    theme_bw() +
    ylab("Effective sample size for historical data") +
    xlab("Bayesian p-value (new vs historical data)")

  op <- par(ask=TRUE)
  plot(post_typeplot)
  plot(densityplot)
  plot(discountfun_plot)
  par(op)
})

#' summary
#'
#' summary
#'
#' @title summary: summary
#' @param object bdpbinomial1arm
#'
#' @examples
#'
#' @method summary
#'
#' @rdname summary
#' @export summary
setMethod("summary", signature(object = "bdpbinomial1arm"), function(object){

  f                   <- object$f1
  posterior_treatment <- object$posterior_treatment

  ### Print
  prior_for_treatment_group <- list(
    `Sample size of prior (treatment group)`           = posterior_treatment$N0,
    `Effective sample size of prior (treatment group)` = posterior_treatment$N0_effective,
    `Bayesian p-value (new vs historical data)`        = posterior_treatment$pvalue,
    `Discount function value`                          = posterior_treatment$alpha_discount)

  ### Text outputs
  print(prior_for_treatment_group)
  argsdf <- data.frame(t(data.frame(object$args1)))
  names(argsdf) <- "args"
  print(argsdf)
})
