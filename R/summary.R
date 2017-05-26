#' @title bdpnormal Object Summary
#' @description \code{summary} method for class \code{bdpnormal}.
#' @param object object of class \code{bdpnormal}. The result of a call to the
#'   \code{\link{bdpnormal}} function.
#'
#' @details Displays a summary of the \code{bdpnormal} fit including the
#'   input data, the stochastic comparison between current and historical
#'   data, and the resulting historical data weight (alpha). If historical
#'   data is missing then no stochastic comparison nor weight are displayed.
#'
#'   In the case of a one-arm analysis, the displayed 95 percent
#'   confidence interval contains the lower and upper limits of the
#'   augmented mean value of the current data. The displayed
#'   \code{mean of treatment group} is the mean of the current data
#'   augmented by the historical data.
#'
#'   When a control arm is present, a two-arm analysis is carried out.
#'   Now, the displayed 95 percent confidence interval contains the
#'   lower and upper limits of the difference between the treatment and
#'   control arms with the historical data augmented to current data, if
#'   present. The displayed augmented sample estimates are the
#'   mean of the treatment and control arms, each of
#'   which are augmented when historical data are present.
#'
#' @import methods
#' @importFrom utils head
#' @importFrom utils write.table
#' @importFrom stats sd density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @export
setMethod("summary", signature(object = "bdpnormal"), function(object){
  arm2                <- object$args$arm2
  posterior_treatment <- object$posterior_treatment
  posterior_control   <- object$posterior_control
  two_side            <- object$args1$two_side
  mu_t                <- object$args1$mu_t
  sigma_t             <- object$args1$sigma_t
  N_t                 <- object$args1$N_t
  mu0_t               <- object$args1$mu0_t
  sigma0_t            <- object$args1$sigma0_t
  N0_t                <- object$args1$N0_t
  mu_c                <- object$args1$mu_c
  sigma_c             <- object$args1$sigma_c
  N_c                 <- object$args1$N_c
  mu0_c               <- object$args1$mu0_c
  sigma0_c            <- object$args1$sigma0_c
  N0_c                <- object$args1$N0_c

  if(!arm2){
    # Format treatment mean output
    mean_CI_t  <- round(quantile(posterior_treatment$posterior_mu, prob=c(0.025, 0.975)),4)
    mean_est_t <- round(median(posterior_treatment$posterior_mu),4)

    cat("\n")
    cat("    One-armed bdp normal\n\n")
    cat("data:\n")
    cat(paste0("  Current treatment: mu_t = ",mu_t,", sigma_t = ", sigma_t, ", N_t = ", N_t))
    cat("\n")
    if(!is.null(N0_t)){
      cat(paste0("  Historical treatment: mu0_t = ",mu0_t,", sigma0_t = ",
                 sigma0_t, ", N0_t = ", N0_t))
      cat("\n")
      cat(paste0("Stochastic comparison (p_hat) - treatment (current vs. historical data): ",
                 round(posterior_treatment$p_hat,4)))
      cat("\n")
      cat(paste0("Discount function value (alpha) - treatment: ",
                 round(posterior_treatment$alpha_discount,4)))
      cat("\n")
    }
    cat("95 percent confidence interval: \n")
    cat(paste0(" ",mean_CI_t))
    cat("\n")
    cat("augmented sample estimate:\n")
    cat("mean of treatment group\n")
    cat(paste0(" ",mean_est_t))
    cat("\n")
  } else{
    comparison_est <- posterior_treatment$posterior_mu - posterior_control$posterior_mu

    mean_est_t <- format(round(median(posterior_treatment$posterior_mu),2), nsmall = 2)
    mean_est_c <- format(round(median(posterior_control$posterior_mu),2), nsmall = 2)
    comp_CI    <- round(quantile(comparison_est, prob=c(0.025, 0.975)),4)
    cat("\n")
    cat("    Two-armed bdp normal\n\n")
    cat("data:\n")
    cat(paste0("  Current treatment: mu_t = ",mu_t,", sigma_t = ", sigma_t, ", N_t = ", N_t))
    cat("\n")
    if(!is.null(N_c)){
      cat(paste0("  Current control: mu_c = ",mu_c,", sigma_c = ", sigma_c, ", N_c = ", N_c))
      cat("\n")
    }
    if(!is.null(N0_t)){
      cat(paste0("  Historical treatment: mu0_t = ",mu0_t,", sigma0_t = ", sigma0_t, ", N0_t = ", N0_t))
      cat("\n")
    }
    if(!is.null(N0_c)){
      cat(paste0("  Historical control: mu0_c = ",mu0_c,", sigma0_c = ", sigma0_c, ", N0_c = ", N0_c))
      cat("\n")
    }

    if(!is.null(N0_t)){
      cat(paste0("Stochastic comparison (p_hat) - treatment (current vs. historical data): ",
                 round(posterior_treatment$p_hat,4)))
      cat("\n")
    }

    if(!is.null(N0_c) & !is.null(N_c)){
      cat(paste0("Stochastic comparison (p_hat) - control (current vs. historical data): ",
                 round(posterior_control$p_hat,4)))
      cat("\n")
    }

    if(!is.null(N0_t)){
      cat(paste0("Discount function value (alpha) - treatment: ",
                 round(posterior_treatment$alpha_discount,4)))
      cat("\n")
    }

    if(!is.null(N0_c) & !is.null(N_c)){
      cat(paste0("Discount function value (alpha) - control: ",
                 round(posterior_control$alpha_discount,4)))
      cat("\n")
    }
    cat("alternative hypothesis: ", ifelse(two_side, "two.sided", "one.sided"))
    cat("\n")

    cat("95 percent confidence interval: \n")
    cat(paste0(" ",comp_CI))
    cat("\n")
    cat("augmented sample estimates:\n")
    cat("treatment group  control group\n")
    cat(paste0("          ", mean_est_t, "          ", mean_est_c))
    cat("\n")
  }

})


#' @title bdpbinomial Object Summary
#' @description \code{summary} method for class \code{bdpbinomial}.
#' @param object object of class \code{bdpbinomial}. The result of a call to the
#'   \code{\link{bdpbinomial}} function.
#'
#' @details Displays a summary of the \code{bdpbinomial} fit including the
#'   input data, the stochastic comparison between current and historical
#'   data, and the resulting historical data weight (alpha). If historical
#'   data is missing then no stochastic comparison nor weight are displayed.
#'
#'   In the case of a one-arm analysis, the displayed 95 percent
#'   confidence interval contains the lower and upper limits of the
#'   augmented event rate of the current data. The displayed
#'   \code{probability of success} is the event rate of the current data
#'   augmented by the historical data.
#'
#'   When a control arm is present, a two-arm analysis is carried out.
#'   Now, the displayed 95 percent confidence interval contains the
#'   lower and upper limits of the difference between the treatment and
#'   control arms with the historical data augmented to current data, if
#'   present. The displayed augmented sample estimates are the
#'   event rates of the treatment and control arms, each of
#'   which are augmented when historical data are present.
#'
#'
#' @import methods
#' @importFrom utils head
#' @importFrom utils write.table
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @export
setMethod("summary", signature(object = "bdpbinomial"), function(object){
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
    if(!is.null(N0_t)){
      cat(paste0("Historical treatment data: ", y0_t, " and ", N0_t, "\n"))
      cat(paste0("Stochastic comparison (p_hat) - treatment (current vs. historical data): ",
                 round(posterior_treatment$p_hat,4)))
      cat("\n")
      cat(paste0("Discount function value (alpha) - treatment: ",
                 round(posterior_treatment$alpha_discount,4)))
      cat("\n")
    }
    cat("95 percent confidence interval: \n")
    cat(paste0(" ",mean_CI_t))
    cat("\n")
    if(!is.null(N0_t)){
      cat("augmented sample estimate:\n")
    } else{
      cat("sample estimates:\n")
    }
    cat("probability of success\n")
    cat(paste0(" ",mean_est_t))
    cat("\n")

  } else{


    ###Format control mean output
    comparison_est <- posterior_treatment$posterior - posterior_control$posterior

    mean_est_t   <- format(round(median(posterior_treatment$posterior),2), nsmall = 2)
    mean_est_c   <- format(round(median(posterior_control$posterior),2), nsmall = 2)

    mean_CI_c    <- round(quantile(posterior_control$posterior, prob=c(0.025, 0.975)),4)
    mean_print_c <- paste0(mean_est_c, " (",mean_CI_c[1], ", ", mean_CI_c[2],")")

    ###Format comparison output
    comp_CI    <- round(quantile(comparison_est, prob=c(0.025, 0.975)),4)
    comp_est   <- round(median(comparison_est),4)
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
      cat(paste0("Stochastic comparison (p_hat) - treatment (current vs. historical data): ",
                 round(posterior_treatment$p_hat,4)))
      cat("\n")
    }
    if(!is.null(N_c) & !is.null(N0_c)){
      cat(paste0("Stochastic comparison (p_hat) - control (current vs. historical data): ",
                 round(posterior_control$p_hat,4)))
      cat("\n")
    }

    if(!is.null(N_t) & !is.null(N0_t)){
      cat(paste0("Discount function value (alpha) - treatment: ",
                 round(posterior_treatment$alpha_discount,4)))
      cat("\n")
    }

    if(!is.null(N_c) & !is.null(N0_c)){
      cat(paste0("Discount function value (alpha) - control: ",
                 round(posterior_control$alpha_discount,4)))
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

  }

})


#' @title bdpsurvival Object Summary
#' @description \code{summary} method for class \code{bdpsurvival}.
#' @import methods
#' @importFrom utils head
#' @importFrom utils write.table
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @param object an object of class \code{bdpsurvival}, a result of a call
#'   to \code{\link{bdpsurvival}}.
#' @details Displays a summary of the \code{bdpsurvival} fit. The output
#'   is different, conditional on a one- or two-arm survival analysis.
#'
#'   In the case of a one-arm analysis, the stochastic comparison between
#'   current and historical data and the resulting historical data weight
#'   (alpha) are displayed, followed by a survival table containing
#'   augmented posterior estimates of the survival probability at each
#'   time point for the current data.
#'
#'   When a control arm is present, a two-arm analysis is carried out.
#'   A two-arm survival analysis is similar to a cox proportional
#'   hazards analysis, and the displayed summary reflects this. First,
#'   the stochastic comparison between current and historical data and
#'   the resulting historical data weight (alpha) are displayed, when
#'   historical data is present for the respective arm. The displayed
#'   \code{coef} value is the log-hazard between the augmented treatment
#'   and control arms (log(treatment) - log(control)). The lower and upper
#'   95 percent interval limits correspond to the \code{coef} value.
#'
#' @export
setMethod("summary", signature(object = "bdpsurvival"), function(object){
  posterior_treatment <- object$posterior_treatment
  posterior_control   <- object$posterior_control
  surv_time           <- object$args1$surv_time

  args1               <- object$args1
  data                <- args1$data
  breaks              <- args1$breaks
  arm2                <- args1$arm2

  historical <- NULL
  treatment  <- NULL

  ##############################################################################
  ### Survival table
  ### - Only print if !arm2
  ##############################################################################
  if(!arm2){
    ### Organize current treatment posterior data
    data_t <- subset(data, historical==0 & treatment == 1)
    s_t    <- with(data_t, Surv(time, status))# , type="mstate"))
    n      <- nrow(data_t)
    s_t    <- survival::survfitKM(factor(rep(1,n)), s_t)

    survival_times_posterior  <- lapply(s_t$time, ppexp,
      posterior_treatment$posterior_hazard,cuts=c(0,breaks))

    s_t$surv    <- 1-sapply(survival_times_posterior, median)
    s_t$std.err <- sapply(survival_times_posterior, sd)
    s_t$upper   <- 1-sapply(survival_times_posterior, quantile, 0.025)
    s_t$lower   <- 1-sapply(survival_times_posterior, quantile, 0.975)

    m_t <- round(cbind(s_t$time, s_t$n.risk, s_t$n.event, s_t$surv, s_t$std.err,
                 s_t$lower, s_t$upper),4)
    cnames <- c("time","n.risk","n.event","survival","std.err",
                "lower 95% CI", "upper 95% CI")
    dimnames(m_t) <- list(rep("", nrow(m_t)), cnames)

    cat("\n")
    cat("    One-armed bdp survival\n\n")
    if(is.null(args1$S0_t)){
      cat("  Current treatment summary:")
      cat("\n")
      print(m_t)
    } else{
      cat(paste0("Stochastic comparison (p_hat) - treatment (current vs. historical data): ",
                 round(posterior_treatment$p_hat,4)))
      cat("\n")
      cat(paste0("Discount function value (alpha) - treatment: ",
                 round(posterior_treatment$alpha_discount,4)))
      cat("\n")
      cat("\n")
      cat("Current treatment - augmented posterior summary:")
      cat("\n")
      print(m_t)
      cat("\n")
    }
  }

  ##############################################################################
  ### Significance of cox proportional hazard for treatment vs control
  ### - Only print if arm2
  ##############################################################################
  if(arm2){
    ### Compute treatment effect of treatment vs control and create table
    R0      <- log(posterior_treatment$posterior_hazard)-log(posterior_control$posterior_hazard)
    V0      <- 1/apply(R0,2,var)
    logHR0  <- R0%*%V0/sum(V0)
    coef    <- mean(logHR0)
    se_coef <- sd(logHR0)
    CI      <- quantile(logHR0, probs=c(0.025,0.975))

    summ_table <- matrix(round(c(coef, exp(coef), se_coef, CI[1], CI[2]),4),nrow=1)
    cnames     <- c("coef", "exp(coef)", "se(coef)", "lower 95% CI", "upper 95% CI")
    dimnames(summ_table) <- list("treatment", cnames)


    ### Count sample size and number of events
    data_t <- subset(data, historical == 0 & treatment == 1)
    n_t    <- nrow(data_t)
    e_t    <- sum(data_t$status)

    if(!is.null(args1$S_c)){
      data_c <- subset(data, historical == 0 & treatment == 0)
      n_c    <- nrow(data_c)
      e_c    <- sum(data_c$status)
    } else{
      data_c <- subset(data, historical == 1 & treatment == 0)
      n_c    <- nrow(data_c)
      e_c    <- sum(data_c$status)
    }

    cat("\n")
    cat("    Two-armed bdp survival\n\n")
    cat("data:\n")
    cat(paste0("  Current treatment: n = ",n_t,", number of events = ", e_t))
    cat("\n")
    if(!is.null(args1$S_c)){
      cat(paste0("  Current control: n = ",n_c,", number of events = ", e_c))
      cat("\n")
    } else{
      cat(paste0("  Historical control: n = ",n_c,", number of events = ", e_c))
      cat("\n")
    }

    if(!is.null(args1$S0_t)){
      cat(paste0("Stochastic comparison (p_hat) - treatment (current vs. historical data): ",
                 round(posterior_treatment$p_hat,4)))
      cat("\n")
    }

    if(!is.null(args1$S0_c) & !is.null(args1$S_c)){
      cat(paste0("Stochastic comparison (p_hat) - control (current vs. historical data): ",
                 round(posterior_control$p_hat,4)))
      cat("\n")
    }

    if(!is.null(args1$S0_t)){
      cat(paste0("Discount function value (alpha) - treatment: ",
                 round(posterior_treatment$alpha_discount,4)))
      cat("\n")
    }

    if(!is.null(args1$S0_c) & !is.null(args1$S_c)){
      cat(paste0("Discount function value (alpha) - control: ",
                 round(posterior_control$alpha_discount,4)))
      cat("\n")
    }
    cat("\n")
    print(summ_table)

  }
})

# Helper functions:

pp <- function(m){
  write.table(format(m, justify="right"),
              row.names=T, col.names=F, quote=F)
}


#' @title bdpregression Object Summary
#' @description \code{summary} method for class \code{bdpregression}.
#' @param object object of class \code{bdpregression}. The result of a call to the
#'   \code{\link{bdpregression}} function.
#'
#' @details Displays a summary of the \code{bdpregression} fit including the
#'   input data, the stochastic comparison between current and historical
#'   data, and the resulting historical data weight (alpha). If historical
#'   data is missing then no stochastic comparison nor weight are displayed.
#'
#'   ...
#'
#' @import methods
#' @importFrom utils head
#' @importFrom utils write.table
#' @export
setMethod("summary", signature(object = "bdpregression"), function(object){

  #arm2                <- object$args$arm2
  #if(!arm2){
    summary(object$posterior_treatment$posterior_regression)
  #s}

})

