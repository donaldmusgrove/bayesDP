#' print
#' print
#' @title print: print
#' @param x bdpnormal
#' @rdname print
#' @aliases print,bdpnormal-method
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


#' print
#' print
#' @title print: print
#' @param x bdpbinomial
#' @rdname print
#' @aliases print,bdpbinomial-method
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


#' print
#' print
#' @title print: print
#' @param x bdpregression_linear
#' @rdname print
#' @aliases print,bdpregression_linear-method
#' @export print
setMethod("print", signature(x = "bdpregression_linear"), function(x){
  f          <- x$f1
  posterior  <- x$est

  ### Print
  prior <- list(`Bayesian p-value (new vs historical data)`       = posterior$pvalue,
                `Loss function value`                             = posterior$alpha_loss)

  print(cat(hypothesis))
  print(prior)
})


#' print
#' print
#' @title print: print
#' @param x bdpsurvival
#' @rdname print
#' @aliases print,bdpsurvival-method
#' @export print
setMethod("print", signature(x = "bdpsurvival"), function(x){
  f                   <- x$f1
  posterior_treatment <- x$posterior_treatment
  surv_time           <- x$args1$surv_time

  ### Print
  prior_for_treatment_group <- list(`Bayesian p-value (new vs historical data)`       = posterior_treatment$pvalue,
                                    `Loss function value`                             = posterior_treatment$alpha_discount)


  ### Output survival probability at requested time
  surv_prob <- list(`Survival time` = surv_time,
                    `Median survival probability` = median(f$treatmentpost))

  ### Text outputs
  print(prior_for_treatment_group)
  print(surv_prob)
})
