#' This code is used for estimating posterior samples from a binary outcome
#' where an informative prior is used. The prior weight is determined using a
#' discount function. In addition this code simulate many trials in order to get
#' trial characteristics you must specify the parameters of the discount function
#' as well as the maximum strength for the prior. This code assumes a
#' non-adaptive trial.
#' This code is modeled after the methodologies developed by the MDIC working
#' group: "Informing clinical trials using bench & simulations"
#' Section 1: of the code defines functions needed
#' Section 2: of the code estimates a posterior and the discount function value
#'            given inputs
#' Section 3: of the code simulates the trial many times to get trial
#'            characteristics
#' Developer: Tarek Haddad
#' Tarek.D.Haddad@Medtronic.com
#' Last modified:1/26/2016
#'
#' bdpbinomial2arm
#'
#' @title bdpbinomial2arm: bdpbinomial2arm
#' @param y_t numeric
#' @param N_t numeric
#' @param y_c numeric
#' @param N_c numeric
#' @param y0_t numeric
#' @param N0_t numeric
#' @param y0_c numeric
#' @param N0_c numeric
#' @param alpha_max numeric
#' @param a0 numeric
#' @param b0 numeric
#' @param number_mcmc numeric
#' @param weibull_shape numeric
#' @param weibull_scale numeric
#' @param two_side character
#' @param inequality character
#' @param delta character
#'
#' @examples
#'
#' @rdname bdpbinomial2arm
#' @export bdpbinomial2arm

setGeneric("bdpbinomial2arm",
           function(y_t           = NULL,
                    N_t           = NULL,
                    y_c           = NULL,
                    N_c           = NULL,
                    y0_t          = NULL,
                    N0_t          = NULL,
                    y0_c          = NULL,
                    N0_c          = NULL,
                    alpha_max     = 1,
                    a0            = 1,
                    b0            = 1,
                    number_mcmc   = 10000,
                    weibull_scale = 0.135,
                    weibull_shape = 3,
                    two_side      = 1){
             standardGeneric("bdpbinomial2arm")
           })

setMethod("bdpbinomial2arm",
          signature(),
          function(y_t           = NULL,
                   N_t           = NULL,
                   y_c           = NULL,
                   N_c           = NULL,
                   y0_t          = NULL,
                   N0_t          = NULL,
                   y0_c          = NULL,
                   N0_c          = NULL,
                   alpha_max     = 1,
                   a0            = 1,
                   b0            = 1,
                   number_mcmc   = 10000,
                   weibull_scale = 0.135,
                   weibull_shape = 3,
                   two_side      = 1){       #Difference margin


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

  f1 <- final_binomial(posterior_treatment = posterior_treatment,
                       posterior_control   = posterior_control)


  args1 <- list(y_t           = y_t,
                N_t           = N_t,
                y_c           = y_c,
                N_c           = N_c,
                y0_t          = y0_t,
                N0_t          = N0_t,
                y0_c          = y0_c,
                N0_c          = N0_c,
                alpha_max     = alpha_max,
                a0            = a0,
                b0            = b0,
                number_mcmc   = number_mcmc,
                weibull_scale = weibull_scale,
                weibull_shape = weibull_shape,
                two_side      = two_side,
                inequality    = inequality,
                delta         = delta)

me <- list(posterior_treatment = posterior_treatment,
           posterior_control   = posterior_control,
           f1                  = f1,
           args1               = args1)

class(me) <- "bdpbinomial2arm"

return(me)

})


#' plot
#'
#' plot
#'
#' @title plot: plot
#' @param x bdpbinomial2arm
#'
#' @examples
#'
#' @rdname plot
#' @export plot
setMethod("plot", signature(x = "bdpbinomial2arm"), function(x){

  f                   <- x$f1
  posterior_treatment <- x$posterior_treatment
  posterior_control   <- x$posterior_control
  two_side            <- x$args1$two_side
  N0_t                <- x$args1$N0_t
  N0_c                <- x$args1$N0_c

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
#' @param x bdpbinomial2arm
#'
#' @examples
#'
#' @rdname print
#' @export print
setMethod("print", signature(x = "bdpbinomial2arm"), function(x){

  f                   <- x$f1
  posterior_treatment <- x$posterior_treatment
  posterior_control   <- x$posterior_control
  two_side            <- x$args1$two_side
  N0_t                <- x$args1$N0_t
  N0_c                <- x$args1$N0_c

  
  if(posterior_treatment$N0==0){
    prior_for_treatment_group <- "No Prior Supplied"
  } else{
    prior_for_treatment_group <- list(
      `Sample size of prior (for treatment group)`          = posterior_treatment$N0,
      `Effective sample size of prior(for treatment group)` = posterior_treatment$N0_effective,
      `Bayesian p-value (new vs historical data)`           = posterior_treatment$pvalue,
      `Discount function value`                             = posterior_treatmentt$alpha_discount)
  }

  if(posterior_control$N0==0){
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
})

#' summary
#'
#' summary
#'
#' @title summary: summary
#' @param object bdpbinomial2arm
#'
#' @examples
#'
#' @rdname summary
#' @export summary
setMethod("summary", signature(object = "bdpbinomial2arm"), function(object){

  f                   <- x$f1
  posterior_treatment <- x$posterior_treatment
  posterior_control   <- x$posterior_control
  two_side            <- x$args1$two_side
  N0_t                <- x$args1$N0_t
  N0_c                <- x$args1$N0_c

  
  if(posterior_treatment$N0==0){
    prior_for_treatment_group <- "No Prior Supplied"
  } else{
    prior_for_treatment_group <- list(
      `Sample size of prior (for treatment group)`          = posterior_treatment$N0,
      `Effective sample size of prior(for treatment group)` = posterior_treatment$N0_effective,
      `Bayesian p-value (new vs historical data)`           = posterior_treatment$pvalue,
      `Discount function value`                             = posterior_treatmentt$alpha_discount)
  }

  if(posterior_control$N0==0){
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
  print(argsdf)
})
