#' summary
#' summary
#' @title summary: summary
#' @importFrom utils head
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @param object Result
#' @export
#' @rdname summary-methods
#' @aliases summary,bdpnormal,bdpnormal-method
setMethod("summary", signature(object = "bdpnormal"), function(object){
  f <- object$f1
  posterior_treatment <- object$posterior_treatment
  posterior_control <- object$posterior_control
  two_side <- object$args1$two_side
  N0_t <- object$args1$N0_t
  N0_c <- object$args1$N0_c

  if(is.null(N0_t) == FALSE){

      prior_for_treatment_group <- list("Sample size of prior (for treatment group)" = N0_t,
                                        "Bayesian p-value (new vs historical data)"  = posterior_treatment$pvalue,
                                        "Discount function value"                    = posterior_treatment$alpha_discount)
      }

  if(is.null(N0_c) == FALSE){
      prior_for_control_group <- list("Sample size of prior (for control group)"  = N0_c,
                                      "Bayesian p-value (new vs historical data)" = posterior_control$pvalue,
                                      "Discount function value"                   = posterior_control$alpha_discount)
      }

  pp(prior_for_treatment_group)
  if(is.null(N0_c) == FALSE){
    pp(prior_for_control_group)
  }
})


#' summary
#' summary
#' @title summary: summary
#' @importFrom utils head
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @export
#' @rdname summary-methods
#' @aliases summary,bdpbinomial,bdpbinomial-method
setMethod("summary", signature(object = "bdpbinomial"), function(object){
  f                   <- object$f1
  posterior_treatment <- object$posterior_treatment
  posterior_control   <- object$posterior_control
  two_side            <- object$args1$two_side
  N0_t                <- object$args1$N0_t
  N0_c                <- object$args1$N0_c

  if(is.null(N0_t) == FALSE){
    prior_for_treatment_group <- list(
      `Sample size of prior (for treatment group)`          = posterior_treatment$N0,
      `Effective sample size of prior(for treatment group)` = posterior_treatment$N0_effective,
      `Bayesian p-value (new vs historical data)`           = posterior_treatment$pvalue,
      `Discount function value`                             = posterior_treatment$alpha_discount)
  }

  if(is.null(N0_t) == FALSE){
    prior_for_control_group <- list(
      `Sample size of prior (for control group)`          = posterior_control$N0,
      `Effective sample size of prior(for control group)` = posterior_control$N0_effective,
      `Bayesian p-value (new vs historical data)`         = posterior_control$pvalue,
      `Discount function value`                           = posterior_control$alpha_discount)
  }

  pp(prior_for_treatment_group)
  if(is.null(posterior_control$N0) == FALSE){
    pp(prior_for_control_group)
  }
})


#' summary
#' summary
#' @title summary: summary
#' @importFrom utils head
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @export
#' @rdname summary-methods
#' @aliases summary,bdpsurvival,bdpsurvival-method
setMethod("summary", signature(object = "bdpsurvival"), function(object){
  f                   <- object$f1
  posterior_treatment <- object$posterior_treatment
  surv_time           <- object$args1$surv_time

  ### Output list
  prior_for_treatment_group <- list("Bayesian p-value (new vs historical data)" = posterior_treatment$pvalue,
                                    "Discount function value"                   = posterior_treatment$alpha_discount,
                                    "Survival time"                             = surv_time,
                                    "Median survival probability"               = median(f$treatmentpost))

  ### Text outputs
  return(prior_for_treatment_group)
})

# Helper functions:

pp <- function(m){
  write.table(format(m, justify="right"),
              row.names=T, col.names=F, quote=F)
}

