#' @title print
#' @description normal print method
#' @import methods
#' @importFrom utils head
#' @importFrom utils write.table
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @param x Result
#' @export
setMethod("print", signature(x = "bdpnormal"), function(x){
  posterior_treatment <- x$posterior_treatment
  posterior_control   <- x$posterior_control
  N0_t                <- x$args1$N0_t
  N0_c                <- x$args1$N0_c

  if(is.null(N0_t) == FALSE){
      prior_for_treatment_group <- list("Sample size of prior (for treatment group):  " = N0_t,
                                        "Bayesian p-value (new vs historical data):  "  = posterior_treatment$p_hat,
                                        "Discount function value:  "                    = posterior_treatment$alpha_discount)
      }
  if(is.null(N0_c) == FALSE){
      prior_for_control_group <- list("Sample size of prior (for control group):  "  = N0_c,
                                      "Bayesian p-value (new vs historical data):  " = posterior_control$p_hat,
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
  invisible(x)
})


#' @title print
#' @description binomial print method
#' @import methods
#' @importFrom utils head
#' @importFrom utils write.table
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @param x Result
#' @export
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
        #`Effective sample size of prior(for treatment group):  ` = posterior_treatment$N0_effective,
        `Bayesian p-value (new vs historical data):  `           = posterior_treatment$p_hat,
        `Discount function value:  `                             = posterior_treatment$alpha_discount)
      }

  if(is.null(N0_c) == FALSE){

      prior_for_control_group <- list(
        `Sample size of prior (for control group):  `          = posterior_control$N0,
        #`Effective sample size of prior(for control group):  ` = posterior_control$N0_effective,
        `Bayesian p-value (new vs historical data):  `         = posterior_control$p_hat,
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
  invisible(x)
})


#' @title print
#' @description survival print method
#' @import methods
#' @importFrom utils head
#' @importFrom utils write.table
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @param x Result
#' @export
setMethod("print", signature(x = "bdpsurvival"), function(x){
  f                   <- x$f1
  posterior_treatment <- x$posterior_treatment
  surv_time           <- x$args1$surv_time
  args1               <- x$args1
  data                <- args1$data
  breaks              <- args1$breaks
  arm2                <- args1$arm2


  if(!arm2){
    ##############################################################################
    # Survival probability and surv_time
    ##############################################################################
    ### Print the augmented posterior
    survival_time_posterior_flat <- ppexp(surv_time,
                                          posterior_treatment$posterior_hazard,
                                          cuts=c(0,breaks))

    data_t <- subset(data, historical==0 & treatment == 1)
    n      <- nrow(data_t)
    s_t    <- with(data_t, Surv(time, status))# , type="mstate"))
    s_t    <- survival::survfitKM(factor(rep(1,n)), s_t)

    print_1arm <- matrix(c(nrow(data_t),
                         sum(s_t$n.event),
                         surv_time,
                         1-median(survival_time_posterior_flat),
                         1-quantile(survival_time_posterior_flat,0.975),
                         1-quantile(survival_time_posterior_flat,0.025)),nrow=1)
    print_1arm <- round(print_1arm,4)
    cnames <- c("n","events","surv_time","median","lower 95% CI","upper 95% CI")
    dimnames(print_1arm) <- list(rep("", nrow(print_1arm)), cnames)
    print(print_1arm)
  }




  # ### Format surv_time output
  # surv_CI    <- round(quantile(f$treatment_posterior, prob=c(0.025, 0.975)),4)
  # surv_est   <- round(median(f$treatment_posterior),4)
  # surv_print <- paste0(surv_est, " (",surv_CI[1], ", ", surv_CI[2],")")
  #
  # ### Output list
  # prior_for_treatment_group <- list(`Stochastic comparison (new vs historical data):  ` = posterior_treatment$p_hat,
  #                                   `Discount function value:  `                        = posterior_treatment$alpha_discount,
  #                                   `Survival time:  `                                  = surv_time,
  #                                   `Median survival probability (95% CI):  `           = surv_print)
  #
  # ### Text outputs
  # cat(pp(prior_for_treatment_group))
  # invisible(x)
})


# Helper functions:

pp <- function(m){
  write.table(format(m, justify="right"),
              row.names=T, col.names=F, quote=F)
}
