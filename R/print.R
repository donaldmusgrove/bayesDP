#' print
#' print
#' @title print: print
#' @importFrom utils head
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @param x Result
#' @export
#' @rdname print-methods
#' @aliases print,bdpnormal,bdpnormal-method
setMethod("print", signature(x = "bdpnormal"), function(x){
  posterior_treatment <- x$posterior_treatment
  posterior_control   <- x$posterior_control
  N0_t                <- x$args1$N0_t
  N0_c                <- x$args1$N0_c

  if(is.null(N0_t) == FALSE){
      prior_for_treatment_group <- list("Sample size of prior (for treatment group):  " = N0_t,
                                        "Bayesian p-value (new vs historical data):  "  = posterior_treatment$pvalue,
                                        "Discount function value:  "                    = posterior_treatment$alpha_discount)
      }
  if(is.null(N0_c) == FALSE){
      prior_for_control_group <- list("Sample size of prior (for control group):  "  = N0_c,
                                      "Bayesian p-value (new vs historical data):  " = posterior_control$pvalue,
                                      "Discount function value:  "                   = posterior_control$alpha_discount)
      }
  pp(prior_for_treatment_group)
  if(is.null(N0_c) == FALSE){
    pp(prior_for_control_group)
  }

  argsdf <- suppressWarnings(data.frame(as.numeric(as.character(x$args1))))
  rownames(argsdf) <- names(x$args1)
  colnames(argsdf) <- "args"
  #print(format(head(argsdf, -2), scientific = FALSE))
  cat("\n")
  cat("Submitted:")
  cat(x$args1$intent)
  cat("\n")
})


#' print
#' print
#' @title print: print
#' @importFrom utils head
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @export
#' @rdname print-methods
#' @aliases print,bdpbinomial,bdpbinomial-method
setMethod("print", signature(x = "bdpbinomial"), function(x){
  f                   <- x$f1
  posterior_treatment <- x$posterior_treatment
  posterior_control   <- x$posterior_control
  two_side            <- x$args1$two_side
  N0_t                <- x$args1$N0_t
  N0_c                <- x$args1$N0_c


  if(is.null(N0_t) == FALSE){
      prior_for_treatment_group <- list(
        `Sample size of prior (for treatment group):  `          = posterior_treatment$N0,
        `Effective sample size of prior(for treatment group):  ` = posterior_treatment$N0_effective,
        `Bayesian p-value (new vs historical data):  `           = posterior_treatment$pvalue,
        `Discount function value:  `                             = posterior_treatment$alpha_discount)
      }

  if(is.null(N0_c) == FALSE){

      prior_for_control_group <- list(
        `Sample size of prior (for control group):  `          = posterior_control$N0,
        `Effective sample size of prior(for control group):  ` = posterior_control$N0_effective,
        `Bayesian p-value (new vs historical data):  `         = posterior_control$pvalue,
        `Discount function value:  `                           = posterior_control$alpha_discount)
      }

  pp(prior_for_treatment_group)
  if(is.null(N0_c) == FALSE){
    pp(prior_for_control_group)
  }

  argsdf <- suppressWarnings(data.frame(as.numeric(as.character(x$args1))))
  rownames(argsdf) <- names(x$args1)
  colnames(argsdf) <- "args"
  print(format(head(argsdf, -2), scientific = FALSE))
  cat("\n")
  cat("Submitted:")
  cat(x$args1$intent)
  cat("\n")
})


#' print
#' print
#' @title print: print
#' @importFrom utils head
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @export
#' @rdname print-methods
#' @aliases print,bdpsurvival,bdpsurvival-method
setMethod("print", signature(x = "bdpsurvival"), function(x){
  f                   <- object$f1
  posterior_treatment <- object$posterior_treatment
  surv_time           <- object$args1$surv_time


  ### Format surv_time output
  surv_CI    <- round(quantile(f$treatmentpost, prob=c(0.025, 0.975)),4)
  surv_est   <- round(median(f$treatmentpost),4)
  surv_print <- paste0(surv_est, " (",surv_CI[1], ", ", surv_CI[2],")")

  ### Output list
  prior_for_treatment_group <- list(`Stochastic comparison (new vs historical data):  ` = posterior_treatment$pvalue,
                                    `Discount function value:  `                        = posterior_treatment$alpha_discount,
                                    `Survival time:  `                                  = surv_time,
                                    `Median survival probability (95% CI):  `           = surv_print)

  ### Text outputs
  return(pp(prior_for_treatment_group))
})


# Helper functions:

pp <- function(m){
  write.table(format(m, justify="right"),
              row.names=T, col.names=F, quote=F)
}
