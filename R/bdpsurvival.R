#' @title Bayesian Discount Prior: Survival Analysis
#' @description \code{bdpsurvival} is used to estimate the survival probability
#'   (single arm trial; OPC) or hazard ratio (two-arm trial; RCT) for
#'   right-censored data using the survival analysis implementation of the
#'   Bayesian discount prior. This code is modeled after
#'   the methodologies developed in Haddad et al. (2017).
#' @param formula an object of class "formula." Must have a survival object on
#'   the left side and exactly two inputs on the right side: treatment and
#'   historical. See "Details" for more information.
#' @param data data frame. A data frame with columns 'time', 'status',
#'   'treatment', and historical.' See "Details" for required structure.
#' @param breaks vector. Breaks (interval starts) used to compose the breaks of the
#'   piecewise exponential model. Do not include zero. Default breaks are the
#'   quantiles of the input times.
#' @param a0 scalar. Prior value for the gamma shape of the piecewise
#'   exponential hazards. Default is 0.1.
#' @param b0 scalar. Prior value for the gamma rate of the piecewise
#'   exponential hazards. Default is 0.1.
#' @param surv_time scalar. Survival time of interest for computing the
#'   probability of survival for a single arm (OPC) trial. Default is
#'   overall, i.e., current+historical, median survival time.
#' @param alpha_max scalar. Maximum weight the discount function can apply.
#'   Default is 1. For a two-arm trial, users may specify a vector of two values
#'   where the first value is used to weight the historical treatment group and
#'   the second value is used to weight the historical control group.
#' @param fix_alpha logical. Fix alpha at alpha_max? Default value is FALSE.
#' @param number_mcmc scalar. Number of Monte Carlo simulations. Default is 10000.
#' @param weibull_shape scalar. Shape parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 3. For a two-arm trial, users may specify a vector of two values
#'   where the first value is used to estimate the weight of the historical
#'   treatment group and the second value is used to estimate the weight of the
#'   historical control group.
#' @param weibull_scale scalar. Scale parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 0.135. For a two-arm trial, users may specify a vector of two
#'   values where the first value is used to estimate the weight of the
#'   historical treatment group and the second value is used to estimate the
#'   weight of the historical control group.
#' @param two_side logical. Indicator of two-sided test for the discount
#'   function. Default value is TRUE.
#' @param method character. Analysis method with respect to estimation of the weight
#'   paramter alpha. Default value "\code{fixed}" estimates alpha once and holds it fixed
#'   throughout the analysis. Alternative method "\code{mc}" estimates alpha for each
#'   Monte Carlo iteration. See the the \code{bdpsurvival} vignette \cr
#'   \code{vignette("bdpsurvival-vignette", package="bayesDP")} for more details.
#' @details \code{bdpsurvival} uses a two-stage approach for determining the
#'   strength of historical data in estimation of a survival probability outcome.
#'   In the first stage, a Weibull distribution function is used as a
#'   \emph{discount function} that defines the maximum strength of the
#'   historical data (via \code{weibull_shape}, \code{weibull_scale}, and
#'   \code{alpha_max}) and discounts based on disagreement with the current data.
#'   Disagreement between current and historical data is determined by stochastically
#'   comparing the respective posterior distributions under noninformative priors.
#'   With a single arm survival data analysis, the comparison is the
#'   probability (\code{p}) that the current survival is less than the historical
#'   survival. For a two-arm survival data, analysis the comparison is the
#'   probability that the hazard ratio comparing treatment and control is
#'   different from zero. The comparison metric \code{p} is then
#'   input into the Weibull discount function and the final strength of the
#'   historical data is returned (alpha).
#'
#'   In the second stage, posterior estimation is performed where the discount
#'   function parameter, \code{alpha}, is used as a fixed value for all posterior
#'   estimation procedures.
#'
#'   To carry out a single arm (OPC) analysis, data for the current and
#'   historical treatments are specified in a data frame. The data frame must have
#'   columns with names 'time', 'status', 'treatment', and 'historical.' Column
#'   'time' is the survival (censor) time of the event and 'status' is the
#'   event indicator. The column 'treatment' is used to indicate which observations
#'   are in the treatment and control group. A value of 1 indicates that the
#'   observation is in the treatment group. The column 'historical' indicates
#'   whether the observation is from the historical data (1) or current data (0).
#'   The results are then based on the posterior probability of survival at
#'   \code{surv_time} for the current data augmented by the historical data.
#'
#'   Two-arm (RCT) analyses are specified similarly to a single arm trial. Again
#'   the input data frame must have columns with names 'time', 'status',
#'   'treatment', and 'historical.' Column 'time' is the survival (censor) time
#'   of the event and 'status' is the event indicator. Now, the 'treatment'
#'   column must use 0 to indicate the control group. The current data are
#'   augmented by historical data (if present) and the results are then based
#'   on the posterior distribution of the hazard ratio between the treatment
#'   and control groups.
#'
#'   For more details, see the \code{bdpsurvival} vignette: \cr
#'   \code{vignette("bdpsurvival-vignette", package="bayesDP")}
#'
#' @return \code{bdpsurvival} returns an object of class "bdpsurvival".
#' The functions \code{\link[=summary,bdpsurvival-method]{summary}} and \code{\link[=print,bdpsurvival-method]{print}} are used to obtain and
#' print a summary of the results, including user inputs. The \code{\link[=plot,bdpsurvival-method]{plot}}
#' function displays visual outputs as well.
#'
#' An object of class "\code{bdpsurvival}" is a list containing at least
#' the following components:
#' \describe{
#'  \item{\code{posterior_treatment}}{
#'    list. Entries contain values related to the treatment group:}
#'    \itemize{
#'      \item{\code{alpha_discount}}{
#'        numeric. Alpha value, the weighting parameter of the historical data.}
#'      \item{\code{p_hat}}{
#'        numeric. The posterior probability of the stochastic comparison
#'        between the current and historical data.}
#'      \item{\code{posterior_survival}}{
#'        vector. If one-arm trial, a vector of length \code{number_mcmc}
#'        containing the posterior probability draws of survival at
#'        \code{surv_time}.}
#'      \item{\code{posterior_flat_survival}}{
#'        vector. If one-arm trial, a vector of length \code{number_mcmc}
#'        containing the probability draws of survival at \code{surv_time}
#'        for the current treatment not augmented by historical treatment.}
#'      \item{\code{prior_survival}}{
#'        vector. If one-arm trial, a vector of length \code{number_mcmc}
#'        containing the probability draws of survival at \code{surv_time}
#'        for the historical treatment.}
#'      \item{\code{posterior_hazard}}{
#'        matrix. A matrix with \code{number_mcmc} rows and \code{length(breaks)}
#'        columns containing the posterior draws of the piecewise hazards
#'        for each interval break point.}
#'      \item{\code{posterior_flat_hazard}}{
#'        matrix. A matrix with \code{number_mcmc} rows and \code{length(breaks)}
#'        columns containing the draws of piecewise hazards for each interval
#'        break point for current treatment not augmented by historical treatment.}
#'      \item{\code{prior_hazard}}{
#'        matrix. A matrix with \code{number_mcmc} rows and \code{length(breaks)}
#'        columns containing the draws of piecewise hazards for each interval break point
#'        for historical treatment.}
#'   }
#'  \item{\code{posterior_control}}{
#'    list. If two-arm trial, contains values related to the control group
#'    analagous to the \code{posterior_treatment} output.}
#'  \item{\code{args1}}{
#'    list. Entries contain user inputs. In addition, the following elements
#'    are ouput:}
#'    \itemize{
#'      \item{\code{S_t}, \code{S_c}, \code{S0_t}, \code{S0_c}}{
#'        survival objects. Used internally to pass survival data between
#'        functions.}
#'      \item{\code{arm2}}{
#'        logical. Used internally to indicate one-arm or two-arm analysis.}
#'   }
#' }
#'
#' @seealso \code{\link[=summary,bdpsurvival-method]{summary}},
#'   \code{\link[=print,bdpsurvival-method]{print}},
#'   and \code{\link[=plot,bdpsurvival-method]{plot}} for details of each of the
#'   supported methods.
#'
#' @references
#' Haddad, T., Himes, A., Thompson, L., Irony, T., Nair, R. MDIC Computer
#'   Modeling and Simulation working group.(2017) Incorporation of stochastic
#'   engineering models as prior information in Bayesian medical device trials.
#'   \emph{Journal of Biopharmaceutical Statistics}, 1-15.
#'
#' @examples
#' # One-arm trial (OPC) example - survival probability at 5 years
#' # Simulate survival data for a single arm (OPC) trial
#' time   <- c(rexp(50, rate=1/20), rexp(50, rate=1/10))
#' status <- c(rexp(50, rate=1/30), rexp(50, rate=1/30))
#' status <- ifelse(time < status, 1, 0)
#'
#' # Collect data into a data frame
#' example_surv_1arm <- data.frame(status     = status,
#'                                 time       = time,
#'                                 historical = c(rep(1,50),rep(0,50)),
#'                                 treatment  = 1)
#'
#' fit1 <- bdpsurvival(Surv(time, status) ~ historical + treatment,
#'                     data = example_surv_1arm,
#'                     surv_time = 5)
#'
#' print(fit1)
#' \dontrun{
#' plot(fit1)
#' }
#'
#' # Two-arm trial (OPC) example
#' # Simulate survival data for a two-arm trial
#' time   <- c(rexp(50, rate=1/20), # Current treatment
#'             rexp(50, rate=1/10), # Current control
#'             rexp(50, rate=1/30), # Historical treatment
#'             rexp(50, rate=1/5))  # Historical control
#' status <- rexp(200, rate=1/40)
#' status <- ifelse(time < status, 1, 0)
#'
#' # Collect data into a data frame
#' example_surv_2arm <- data.frame(status     = status,
#'                                 time       = time,
#'                                 historical = c(rep(0,100),rep(1,100)),
#'                                 treatment  = c(rep(1,50),rep(0,50),rep(1,50),rep(0,50)))
#'
#' fit2 <- bdpsurvival(Surv(time, status) ~ historical + treatment,
#'                        data = example_surv_2arm)
#'
#' summary(fit2)
#'
#' ### Fix alpha at 1
#' fit2_1 <- bdpsurvival(Surv(time, status) ~ historical + treatment,
#'                       data = example_surv_2arm,
#'                       fix_alpha = TRUE)
#'
#' summary(fit2_1)
#'
#'
#' @rdname bdpsurvival
#' @import methods
#' @importFrom stats sd density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov pchisq
#' @importFrom survival Surv survSplit
#' @aliases bdpsurvival,ANY-method
#' @export bdpsurvival
bdpsurvival <- setClass("bdpsurvival", slots = c(posterior_treatment = "list",
                                                 f1 = "list",
                                                 args1 = "list"))
setGeneric("bdpsurvival",
  function(formula       = formula,
           data          = data,
           breaks        = NULL,
           a0            = 0.1,
           b0            = 0.1,
           surv_time     = NULL,
           alpha_max     = 1,
           fix_alpha     = FALSE,
           number_mcmc   = 10000,
           weibull_scale = 0.135,
           weibull_shape = 3,
           two_side      = TRUE,
           method        = "fixed"){
             standardGeneric("bdpsurvival")
           })

setMethod("bdpsurvival",
  signature(),
  function(formula       = formula,
           data          = data,
           breaks        = NULL,
           a0            = 0.1,
           b0            = 0.1,
           surv_time     = NULL,
           alpha_max     = 1,
           fix_alpha     = FALSE,
           number_mcmc   = 10000,
           weibull_scale = 0.135,
           weibull_shape = 3,
           two_side      = TRUE,
           method        = "fixed"){

  ### Check data frame and ensure it has the correct column names
  namesData <- tolower(names(data))
  namesDiff <- setdiff(c("status", "time", "historical", "treatment"), namesData)
  if(length(namesDiff)>0){
    nDiff <- length(namesDiff)
    if(nDiff == 1){
      errorMsg <- paste0("Column ",
                         namesDiff,
                         " is missing from the input data frame.")
      stop(errorMsg)
    } else if(nDiff>1){
      errorNames <- paste0(namesDiff, collapse = ", ")
      errorMsg <- paste0("Columns are missing from input data frame: ",
                         errorNames)
      stop(errorMsg)
    }
  }


  historical <- NULL
  treatment <- NULL

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
  # Format input data
  ##############################################################################
  ### If no breaks input, create intervals along quantiles
  if(is.null(breaks)){
     breaks <- quantile(data$time,probs=c(0.2,0.4,0.6,0.8))
  }

  ### If zero is present in breaks, remove and give warning
  if(any(breaks==0)){
    breaks <- breaks[!(breaks==0)]
    warning("Breaks vector included 0. The zero value was removed.")
  }



  ### Split the data on the breaks
  dataSplit <- survSplit(formula,
                         cut     = breaks,
                         start   = "start",
                         episode = "interval",
                         data    = data)


  ### Compute exposure time within each interval
  dataSplit$exposure <- dataSplit$time - dataSplit$start

  ### Create new labels for the intervals
  maxTime  <- max(dataSplit$time)
  labels_t <- unique(c(breaks, maxTime))
  dataSplit$interval <- factor(dataSplit$interval,
                               labels = labels_t)

  ### Parse out the historical and current data
  S_t  <- subset(dataSplit, historical==0 & treatment == 1)
  S_c  <- subset(dataSplit, historical==0 & treatment == 0)
  S0_t <- subset(dataSplit, historical==1 & treatment == 1)
  S0_c <- subset(dataSplit, historical==1 & treatment == 0)

  if(nrow(S0_t) == 0) S0_t <- NULL
  if(nrow(S_c) == 0)  S_c  <- NULL
  if(nrow(S0_c) == 0) S0_c <- NULL


  ### Compute arm2, internal indicator of a two-arm trial
  if(is.null(S_c) & is.null(S0_c)){
    arm2 <- FALSE
  } else{
    arm2 <- TRUE
  }

  #if(arm2) stop("Two arm trials are not currently supported.")

  ### If surv_time is null, replace with median time
  if(is.null(surv_time) & !arm2){
    surv_time <- median(data$time)
  }


  ### Check inputs
  if(!arm2){
    if(nrow(S_t) == 0) stop("Current treatment data missing or input incorrectly.")
    if(is.null(S0_t)) warning("Historical treatment data missing or input incorrectly.")
  } else if(arm2){
    if(nrow(S_t) == 0) stop("Current treatment data missing or input incorrectly.")
    if(is.null(S_c)) warning("Current control data missing or input incorrectly.")
    if(is.null(S0_t) & is.null(S0_c)) warning("Historical data input incorrectly.")
  }

  posterior_treatment <- posterior_survival(
    S             = S_t,
    S0            = S0_t,
    surv_time     = surv_time,
    alpha_max     = alpha_max[1],
    fix_alpha     = fix_alpha,
    a0            = a0,
    b0            = b0,
    number_mcmc   = number_mcmc,
    weibull_shape = weibull_shape[1],
    weibull_scale = weibull_scale[1],
    two_side      = two_side,
    breaks        = breaks,
    arm2          = arm2,
    method        = method)

  if(arm2){
    posterior_control <- posterior_survival(
      S             = S_c,
      S0            = S0_c,
      alpha_max     = alpha_max[2],
      fix_alpha     = fix_alpha,
      a0            = a0,
      b0            = b0,
      surv_time     = surv_time,
      number_mcmc   = number_mcmc,
      weibull_shape = weibull_shape[2],
      weibull_scale = weibull_scale[2],
      two_side      = two_side,
      breaks        = breaks,
      arm2          = arm2,
      method        = method)
  } else{
    posterior_control <- NULL
  }


  args1 <- list(S_t           = S_t,
                S_c           = S_c,
                S0_t          = S0_t,
                S0_c          = S0_c,
                alpha_max     = alpha_max,
                fix_alpha     = fix_alpha,
                a0            = a0,
                b0            = b0,
                surv_time     = surv_time,
                number_mcmc   = number_mcmc,
                weibull_scale = weibull_scale,
                weibull_shape = weibull_shape,
                two_side      = two_side,
                method        = method,
                arm2          = arm2,
                breaks        = breaks,
                data          = data)

  me <- list(posterior_treatment = posterior_treatment,
             posterior_control   = posterior_control,
             args1               = args1)

  class(me) <- "bdpsurvival"

  return(me)
})




################################################################################
# Survival posterior estimation
# 1) Estimate the discount function (if current+historical data both present)
# 2) Estimate the posterior of the augmented data
################################################################################
### Combine  loss function and posterior estimation into one function
posterior_survival <- function(S, S0, surv_time, alpha_max, fix_alpha, a0, b0,
                               number_mcmc, weibull_shape, weibull_scale,
                               two_side, breaks, arm2, method){


  ### Extract intervals and count number of intervals
  ### - It should be that S_int equals S0_int
  if(!is.null(S)){
    S_int  <- levels(S$interval)
    nInt   <- length(S_int)
  }

  if(!is.null(S0)){
    S0_int  <- levels(S0$interval)
    nInt    <- length(S0_int)
  }

  interval <- NULL


  ##############################################################################
  # Discount function
  # - Comparison is made only if both S and S0 are present
  ##############################################################################
  # Compute hazards for historical and current data efficiently
  if(!is.null(S) & !is.null(S0)){
    ### Compute posterior of interval hazards
    a_post  <- b_post  <- numeric(nInt)
    a_post0 <- b_post0 <- numeric(nInt)

    if(!is.null(S)){
      posterior_flat_hazard <- prior_hazard <- matrix(NA, number_mcmc, nInt)
    } else{
      posterior_flat_hazard <- prior_hazard <- matrix(NA, number_mcmc, nInt)
    }


    ### Compute posterior values
    for(i in 1:nInt){
      a_post[i] <- a0 + sum(subset(S, interval==S_int[i])$status)
      b_post[i] <- b0 + sum(subset(S, interval==S_int[i])$exposure)

      a_post0[i] <- a0 + sum(subset(S0, interval==S0_int[i])$status)
      b_post0[i] <- b0 + sum(subset(S0, interval==S0_int[i])$exposure)

      ### Interval hazards - add on a very small value to avoid underflow
      posterior_flat_hazard[,i]  <- rgamma(number_mcmc, a_post[i],  b_post[i])+1e-12
      prior_hazard[,i]           <- rgamma(number_mcmc, a_post0[i], b_post0[i])+1e-12
    }
  } else if(!is.null(S) & is.null(S0)){
    ### Compute posterior of interval hazards
    a_post  <- b_post  <- numeric(nInt)

    posterior_flat_hazard <- matrix(NA, number_mcmc, nInt)
    prior_hazard          <- NULL

    ### Compute posterior values
    for(i in 1:nInt){
      a_post[i] <- a0 + sum(subset(S, interval==S_int[i])$status)
      b_post[i] <- b0 + sum(subset(S, interval==S_int[i])$exposure)

      ### Interval hazards - add on a very small value to avoid underflow
      posterior_flat_hazard[,i]  <- rgamma(number_mcmc, a_post[i],  b_post[i])+1e-12
    }
  } else if(is.null(S) & !is.null(S0)) {
    ### Compute posterior of interval hazards
    a_post0 <- b_post0 <- numeric(nInt)

    prior_hazard          <- matrix(NA, number_mcmc, nInt)
    posterior_flat_hazard <- NULL

    ### Compute posterior values
    for(i in 1:nInt){
      a_post0[i] <- a0 + sum(subset(S0, interval==S0_int[i])$status)
      b_post0[i] <- b0 + sum(subset(S0, interval==S0_int[i])$exposure)

      ### Interval hazards - add on a very small value to avoid underflow
      prior_hazard[,i]           <- rgamma(number_mcmc, a_post0[i], b_post0[i])+1e-12
    }
  }


  ### If only one of S or S0 is present, return related hazard and (if !arm2), return survival
  if(!is.null(S) & is.null(S0)){
    posterior_hazard <- posterior_flat_hazard

    if(!arm2){
      posterior_survival <- posterior_flat_survival <- 1 - ppexp(q=surv_time,
                                                                 x=posterior_hazard,
                                                                 cuts = c(0,breaks))
      prior_survival      <- NULL
    } else{
      posterior_survival <- posterior_flat_survival <- prior_survival <- NULL
    }
  } else if(is.null(S) & !is.null(S0)){
    posterior_hazard <- prior_hazard

    if(!arm2){
      posterior_survival <- prior_survival <- 1 - ppexp(q=surv_time,
                                                        x=posterior_hazard,
                                                        cuts = c(0,breaks))
      posterior_flat_survival <- NULL
    } else{
      posterior_survival  <- prior_survival <- posterior_flat_survival <- NULL
    }
  }

  if(!(!is.null(S) & !is.null(S0))){
    return(list(alpha_discount          = NULL,
                p_hat                   = NULL,
                posterior_survival      = posterior_survival,
                posterior_flat_survival = posterior_flat_survival,
                prior_survival          = prior_survival,
                posterior_hazard        = posterior_hazard,
                posterior_flat_hazard   = posterior_flat_hazard,
                prior_hazard            = prior_hazard))
  }


  ### If both S and S0 are present, carry out the comparison and compute alpha
  if(!arm2){
    ### Posterior survival probability
    posterior_flat_survival  <- 1 - ppexp(q=surv_time, x=posterior_flat_hazard, cuts = c(0,breaks))
    prior_survival           <- 1 - ppexp(q=surv_time, x=prior_hazard, cuts = c(0,breaks))

    ### Compute probability that survival is greater for current vs historical
    if(method == "mc"){
      logS  <- log(posterior_flat_survival)
      logS0 <- log(prior_survival)

      ### Variance of log survival, computed via delta method of hazards
      nIntervals <- sum(surv_time > c(0,breaks))
      IntLengths <- c(c(0,breaks)[1:nIntervals], surv_time)
      surv_times <- diff(IntLengths)

      v          <- as.matrix(posterior_flat_hazard[,1:nIntervals]^2) %*% (surv_times^2/a_post[1:nIntervals])
      v0         <- as.matrix(prior_hazard[,1:nIntervals]^2) %*% (surv_times^2/a_post0[1:nIntervals])
      Z          <- (logS - logS0)^2 / (v+v0)
      p_hat      <- pchisq(Z, df=1, lower.tail=FALSE)

    } else if(method == "fixed"){
      p_hat <- mean(posterior_flat_survival > prior_survival)   # higher is better survival
    } else{
      stop("Unrecognized method. Use one of 'fixed' or 'mc'")
    }

    if(fix_alpha){
      alpha_discount <- alpha_max
    } else{
      if (!two_side | method == "mc") {
        alpha_discount <- pweibull(p_hat, shape=weibull_shape, scale=weibull_scale)*alpha_max
      } else if (two_side){
        p_hat    <- ifelse(p_hat > 0.5, 1 - p_hat, p_hat)
        alpha_discount <- pweibull(p_hat, shape=weibull_shape, scale=weibull_scale)*alpha_max
      }
    }
  } else{
    ### Weight historical data via (approximate) hazard ratio comparing
    ### current vs historical
    if(method == "mc"){
      R0    <- log(prior_hazard)-log(posterior_flat_hazard)
      v     <- 1/(a_post)
      v0    <- 1/(a_post0)
      R     <- rowSums(as.matrix(R0/(v+v0) / sum(1/(v+v0))))
      V0    <- 1 / sum(1/(v+v0))
      Z     <- R^2 / V0
      p_hat <- pchisq(Z, df=1, lower.tail=FALSE)
    } else if(method == "fixed"){
      R0     <- log(prior_hazard)-log(posterior_flat_hazard)
      V0     <- 1/apply(R0,2,var)
      logHR0 <- R0%*%V0/sum(V0)    #weighted average  of SE^2
      p_hat <- mean(logHR0 > 0)    #larger is higher failure
    } else{
      stop("Unrecognized method. Use one of 'fixed' or 'mc'")
    }

    if(fix_alpha){
      alpha_discount <- alpha_max
    } else{
      if (!two_side | method == "mc") {
        alpha_discount <- pweibull(p_hat, shape=weibull_shape, scale=weibull_scale)*alpha_max
      } else if (two_side){
        p_hat    <- ifelse(p_hat > 0.5, 1 - p_hat, p_hat)
        alpha_discount <- pweibull(p_hat, shape=weibull_shape, scale=weibull_scale)*alpha_max
      }
    }
  }



  ##############################################################################
  # Posterior augmentation via the interval hazards
  # - If current or historical data are missing, this will not augment(see above)
  ##############################################################################
  posterior_hazard <- matrix(NA, number_mcmc, nInt)

  for(i in 1:nInt){
    a_post0[i] <- a0 + sum(subset(S0, interval==S0_int[i])$status)
    b_post0[i] <- b0 + sum(subset(S0, interval==S0_int[i])$exposure)

    a_post[i] <- a0 + sum(subset(S, interval==S_int[i])$status)
    b_post[i] <- b0 + sum(subset(S, interval==S_int[i])$exposure)

    ### Add on a very small value to avoid underflow
    posterior_hazard[,i]  <- rgamma(number_mcmc,
                                    a_post[i]+a_post0[i]*alpha_discount,
                                    b_post[i]+b_post0[i]*alpha_discount) + 1e-12
  }

  ### Posterior of survival time (if !arm2)
  if(!arm2){
    posterior_survival <- 1-ppexp(q=surv_time, x=posterior_hazard, cuts=c(0,breaks))
  } else{
    posterior_survival <- posterior_flat_survival <- prior_survival <- NULL
  }


  return(list(alpha_discount          = alpha_discount,
              p_hat                   = p_hat,
              posterior_survival      = posterior_survival,
              posterior_flat_survival = posterior_flat_survival,
              prior_survival          = prior_survival,
              posterior_hazard        = posterior_hazard,
              posterior_flat_hazard   = posterior_flat_hazard,
              prior_hazard            = prior_hazard))
}

