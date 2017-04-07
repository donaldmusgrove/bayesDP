#' @title summary
#' @description normal summary method
#' @import methods
#' @importFrom utils head
#' @importFrom utils write.table
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @param object Result
#' @export
#' @rdname summary-methods
#' @aliases summary,bdpnormal
setMethod("summary", signature(object = "bdpnormal"), function(object){
  f <- object$f1
  posterior_treatment <- object$posterior_treatment
  posterior_control <- object$posterior_control
  two_side <- object$args1$two_side
  N0_t <- object$args1$N0_t
  N0_c <- object$args1$N0_c

  if(is.null(N0_t) == FALSE){
    ###Format treatment mean output
    mean_CI_t    <- round(quantile(posterior_treatment$mu_posterior, prob=c(0.025, 0.975)),4)
    mean_est_t   <- round(median(posterior_treatment$mu_posterior),4)
    mean_print_t <- paste0(mean_est_t, " (",mean_CI_t[1], ", ", mean_CI_t[2],")")

    prior_for_treatment_group <- list(
      `Stochastic comparison - treatment (new vs historical data):  ` = posterior_treatment$pvalue,
      `Discount function value - treatment:  `                        = posterior_treatment$alpha_discount,
      `Sample size of historical - treatment:  `                      = N0_t,
      `Mean posterior - treatment (95% CI):  `                        = mean_print_t
    )
  }

  if(is.null(N0_c) == FALSE){
    ###Format control mean output
    mean_CI_c    <- round(quantile(posterior_control$mu_posterior, prob=c(0.025, 0.975)),4)
    mean_est_c   <- round(median(posterior_control$mu_posterior),4)
    mean_print_c <- paste0(mean_est_c, " (",mean_CI_c[1], ", ", mean_CI_c[2],")")

    ###Format comparison output
    comp_CI    <- round(quantile(f$TestMinusControl_post, prob=c(0.025, 0.975)),4)
    comp_est   <- round(median(f$TestMinusControl_post),4)
    comp_print <- paste0(comp_est, " (",comp_CI[1], ", ", comp_CI[2],")")

    prior_for_control_group <- list(
      `Stochastic comparison control (new vs historical data):  ` = posterior_control$pvalue,
      `Discount function value - control:  `                      = posterior_control$alpha_discount,
      `Sample size of historical - control:  `                    = N0_c,
      `Mean posterior - control (95% CI):  `                      = mean_print_c,
      `Mean difference - treatment vs. control (95% CI): `        = comp_print
    )
  }

  pp(prior_for_treatment_group)
  if(is.null(N0_c) == FALSE){
    pp(prior_for_control_group)
  }
})


#' @title summary
#' @description binomial summary method
#' @importFrom utils head
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @export
#' @rdname summary-methods
#' @aliases summary,bdpbinomial
setMethod("summary", signature(object = "bdpbinomial"), function(object){
  f                   <- object$f1
  arm2                <- object$args$arm2
  posterior_treatment <- object$posterior_treatment
  posterior_control   <- object$posterior_control
  two_side            <- object$args1$two_side
  y_t                 <- object$args1$y_t
  N_t                 <- object$args1$N_t
  y_c                 <- object$args1$y_c
  N_c                 <- object$args1$N_c
  y0_t                <- object$args1$y0_t
  N0_t                <- object$args1$N0_t
  y0_c                <- object$args1$y0_c
  N0_c                <- object$args1$N0_c

  if(!arm2){
    ###Format treatment mean output
    mean_CI_t    <- round(quantile(posterior_treatment$posterior, prob=c(0.025, 0.975)),4)
    mean_est_t   <- round(median(posterior_treatment$posterior),4)
    mean_print_t <- paste0(mean_est_t, " (",mean_CI_t[1], ", ", mean_CI_t[2],")")

    cat("\n")
    cat("    One-armed bdp binomial\n\n")
    cat(paste0("Current treatment data: ", y_t, " and ", N_t, "\n"))
    cat(paste0("Historical treatment data: ", y0_t, " and ", N0_t, "\n"))
    cat(paste0("Stochastic comparison - treatment (current vs. historical data): ", posterior_treatment$pvalue))
    cat("\n")
    cat(paste0("Discount function value - treatment: ", posterior_treatment$alpha_discount))
    cat("\n")
    cat("95 percent confidence interval: \n")
    cat(paste0(" ",mean_CI_t))
    cat("\n")
    cat("sample estimates:\n")
    cat("probability of success\n")
    cat(paste0(" ",mean_est_t))
    cat("\n")

    #prior_for_treatment_group <- list(
    #  `Stochastic comparison - treatment (new vs historical data):  ` = posterior_treatment$pvalue,
    #  `Discount function value - treatment:  `                        = posterior_treatment$alpha_discount,
    #  `Sample size of historical - treatment:  `                      = posterior_treatment$N0,
    #  `Binomial posterior - treatment (95% CI):  `                    = mean_print_t)

  } else{


    ###Format control mean output

    mean_est_t   <- format(round(median(posterior_treatment$posterior),2), nsmall = 2)
    mean_est_c   <- format(round(median(posterior_control$posterior),2), nsmall = 2)

    mean_CI_c    <- round(quantile(posterior_control$posterior, prob=c(0.025, 0.975)),4)
    mean_print_c <- paste0(mean_est_c, " (",mean_CI_c[1], ", ", mean_CI_c[2],")")

    ###Format comparison output
    comp_CI    <- round(quantile(f$comparison_posterior, prob=c(0.025, 0.975)),4)
    comp_est   <- round(median(f$comparison_posterior),4)
    comp_print <- paste0(comp_est, " (",comp_CI[1], ", ", comp_CI[2],")")

    cat("\n")
    cat("    Two-armed bdp binomial\n\n")
    cat(paste0("Current treatment data: ", y_t, " and ", N_t, "\n"))
    if(!is.null(N_c)){
      cat(paste0("Current control data: ", y_c, " and ", N_c, "\n"))
    }
    if(!is.null(N0_t)){
      cat(paste0("Historical treatment data: ", y0_t, " and ", N0_t, "\n"))
    }
    if(!is.null(N0_c)){
      cat(paste0("Historical control data: ", y0_c, " and ", N0_c, "\n"))
    }
    if(!is.null(N_t) & !is.null(N0_t)){
      cat(paste0("Stochastic comparison - treatment (current vs. historical data): ", posterior_treatment$pvalue))
      cat("\n")
    }
    if(!is.null(N_c) & !is.null(N0_c)){
      cat(paste0("Stochastic comparison - control (current vs. historical data): ", posterior_control$pvalue))
      cat("\n")
    }

    if(!is.null(N_t) & !is.null(N0_t)){
      cat(paste0("Discount function value - treatment: ", posterior_treatment$alpha_discount))
      cat("\n")
    }

    if(!is.null(N_c) & !is.null(N0_c)){
      cat(paste0("Discount function value - control: ", posterior_control$alpha_discount))
      cat("\n")
    }

    cat("alternative hypothesis: ", ifelse(two_side, "two.sided", "one.sided"))
    cat("\n")
    cat("95 percent confidence interval:\n")
    cat(paste0(" "), comp_CI)
    cat("\n")
    cat("augmented sample estimates:\n")
    cat("prop 1 prop2\n")
    cat(paste0("  ", mean_est_t, "  ", mean_est_c))


    #prior_for_control_group <- list(
    #  `Stochastic comparison - control (new vs historical data):  ` = posterior_control$pvalue,
    #  `Discount function value - control:  `                        = posterior_control$alpha_discount,
    #  `Sample size of historical - control:  `                      = posterior_control$N0,
    #  `Binomial posterior - control (95% CI):  `                    = mean_print_c,
    #  `Binomial difference - treatment vs. control (95% CI): `      = comp_print)
  }

  #pp(prior_for_treatment_group)
  #if(is.null(posterior_control$N0) == FALSE){
  #  pp(prior_for_control_group)
  #}
})


#' @title summary
#' @description survival summary method
#' @importFrom utils head
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @export
#' @rdname summary-methods
#' @aliases summary,bdpsurvival
setMethod("summary", signature(object = "bdpsurvival"), function(object){
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

