#' summary
#' summary
#' @title summary: summary
#' @param object bdpnormal
#' @rdname summary
#' @aliases summary,bdpnormal-method
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
  print(format(head(argsdf, -2), scientific = FALSE))
  cat("\n")
  cat("Submitted:")
  cat(object$args1$intent)
})


#' summary
#' summary
#' @title summary: summary
#' @param object bdpbinomial
#' @rdname summary
#' @aliases summary,bdpbinomial-method
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
  if(is.null(posterior_control$N0) == FALSE){
    print(prior_for_control_group)
  }

  argsdf <- suppressWarnings(data.frame(as.numeric(as.character(object$args1))))
  rownames(argsdf) <- names(object$args1)
  colnames(argsdf) <- "args"
  print(format(head(argsdf, -2), scientific = FALSE))
  cat("\n")
  cat("Submitted:")
  cat(object$args1$intent)
})


#' summary
#' summary
#' @title summary: summary
#' @param object bdpregression_linear
#' @rdname summary
#' @aliases summary,bdpregression_linear-method
#' @export summary
setMethod("summary", signature(object = "bdpregression_linear"), function(object){
  f          <- object$f1
  posterior  <- object$est

  ### summary
  prior <- list(`Bayesian p-value (new vs historical data)`       = posterior$pvalue,
                `Loss function value`                             = posterior$alpha_loss)

  print(cat(hypothesis))
  print(prior)
})


#' summary
#' summary
#' @title summary: summary
#' @param object bdpsurvival
#' @rdname summary
#' @aliases summary,bdpsurvival-method
#' @export summary
setMethod("summary", signature(object = "bdpsurvival"), function(object){
  f                   <- object$f1
  posterior_treatment <- object$posterior_treatment
  surv_time           <- object$args1$surv_time

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
