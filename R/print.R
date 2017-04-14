#' @title print
#' @description normal print method
#' @import methods
#' @importFrom utils head
#' @importFrom utils write.table
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @param x Result
#' @export
setMethod("print", signature(x = "bdpnormal"), function(x){
  ### Return summary
  summary(x)
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
  ### Return summary
  summary(x)
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
  posterior_treatment <- x$posterior_treatment
  posterior_control   <- x$posterior_control
  surv_time           <- x$args1$surv_time

  args1               <- x$args1
  data                <- args1$data
  breaks              <- args1$breaks
  arm2                <- args1$arm2


  treatment = NULL
  historical = NULL

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
    cat("\n")
    cat("    One-armed bdp survival\n\n")
    cat("\n")
    print(print_1arm)
  } else{
    ### Return summary
    summary(x)
  }

})


# Helper functions:

pp <- function(m){
  write.table(format(m, justify="right"),
              row.names=T, col.names=F, quote=F)
}
