#' @title Bayesian Discount Prior: Survival Analysis
#' @description \code{bdpsurvival} is used to estimate the survival probability
#'   (single arm trial; OPC) or hazard ratio (two-arm trial; RCT) for
#'   right-censored data using the survival analysis implementation of the
#'   Bayesian discount prior.
#' @param formula an object of class "formula." Must have a survival object on
#'   the left side and exactly two inputs on the right side: treatment and
#'   historical. See "Details" for more information.
#' @param data data frame. A data frame with columns 'time', 'status',
#'   'treatment', and historical.' See "Details" for required structure.
#' @param breaks vector. Breaks (interval starts) used to compose the breaks of the
#'   piecewise exponential model. Do not include zero.
#' @param a0 scalar. Prior value for the gamma shape. Default is 1.
#' @param b0 scalar. Prior value for the gamma rate. Default is 1.
#' @param surv_time scalar. Survival time of interest for computing the the
#'   probability of survival for a single arm, i.e., an OPC trial. Default is
#'   median survival time.
#' @param alpha_max scalar. Maximum weight the discount function can apply.
#'   Default is 1. For a two-arm trial, users may specify a vector of two values
#'   where the first value is used to weight the historical treatment group and
#'   the second value is used to weight the historical control group.
#' @param fix_alpha logical. Fix alpha at alpha_max? Default value is FALSE.
#' @param number_mcmc scalar. Number of Markov Chain Monte Carlo (MCMC)
#'   simulations. Default is 10000.
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
#'
#' @details \code{bdpsurvival} uses a two-stage approach for determining the
#'   strength of historical data in estimation of a survival probability outcome.
#'   In the first stage, a Weibull distribution function is used as a
#'   \emph{discount function} that defines the maximum strength of the
#'   historical data (via \code{weibull_shape}, \code{weibull_scale}, and
#'   \code{alpha_max}) and discounts based on disagreement with the current data.
#'   Disagreement between current and historical data is determined by stochastically
#'   comparing the respective posterior distributions under noninformative priors.
#'   With a single arm survival data analysis, the comparison is the
#'   proability (\code{p}) that the current survial is less than the historical
#'   survival. The comparison metric \code{p} is then
#'   input into the Weibull discount function and the final strength of the
#'   historical data is returned (alpha).
#'
#' In the second stage, posterior estimation is performed where the discount
#'   function parameter, \code{alpha}, is used as a fixed value for all posterior
#'   estimation procedures.
#'
#'  To carry out a single arm (OPC) analysis, data for the current and
#'  historical treatments are specified in a dataframe. The dataframe must have
#'  columns with names 'time', 'status', 'treatment', and 'historical.' Column
#'  'time' is the survival (censor) time of the event and 'status' is the
#'  event indicator. The column 'treatment' is used to indicate which observations
#'  are in the treatment and control group. A value of 1 indicates that the
#'  observation is in the treatment group. The column 'historical' indicates
#'  whether the observation is from the historical data (1) or current data (0).
#'  The results are then based on the posterior distribution of the current data
#'  augmented by the historical data.
#'
#'  Two-arm (RCT) analyses are not available with this release.
#'
#' @return \code{bdpsurvival} returns an object of class "bdpsurvival".
#' The functions \code{summary} and \code{print} are used to obtain and
#' print a summary of the results, including user inputs. The \code{plot}
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
#'      \item{\code{posterior}}{
#'        list. Entries contain \code{cdf}, a vector of the posterior cdf of the
#'        piecewise exponentialdistribution; \code{surv_time_posterior}, a
#'        vector with the posterior of the survival probability; and
#'        \code{posterior}, a matrix of the posteriors of each of the piecewise
#'        hazards.}
#'      \item{\code{posterior_flat}}{
#'        matrix. The distributions of the current treatment group piecewise
#'        hazard rates.}
#'      \item{\code{prior}}{
#'        matrix. The distributions of the historical treatment group piecewise
#'        hazard rates.}
#'   }
#'  \item{\code{f1}}{
#'    list. Entries contain values related to the posterior effect:}
#'    \itemize{
#'      \item{\code{density_post_treatment}}{
#'        object of class \code{density}. Used internally to plot the density of
#'        the treatment group posterior.}
#'      \item{\code{density_flat_treatment}}{
#'        object of class \code{density}. Used internally to plot the density of
#'        the treatment group "flat" distribution.}
#'      \item{\code{density_prior_treatment}}{
#'        object of class \code{density}. Used internally to plot the density of
#'        the treatment group prior.}
#'      \item{\code{treatmentpost}}{
#'        vector. Used internally to plot the posterior distribution of the
#'        survival probability.}
#'   }
#'  \item{\code{args1}}{
#'    list. Entries contain user inputs. In addition, the following elements
#'    are ouput:}
#'    \itemize{
#'      \item{\code{S_t} and \code{S0_t}}{
#'        survival objects. Used internally pass survival data between
#'        functions.}
#'   }
#' }
#'
#' @examples
#' # One-arm trial (OPC) example
#' # Simulate survival data for a single arm (OPC) trial
#' time   <- c(rexp(50, rate=1/20), rexp(50, rate=1/10))
#' status <- c(rexp(50, rate=1/30), rexp(50, rate=1/30))
#' status <- ifelse(time < status, 1, 0)
#'
#' # Collect data into a dataframe
#' example_surv_1arm <- data.frame(status     = status,
#'                                 time       = time,
#'                                 historical = c(rep(1,50),rep(0,50)),
#'                                 treatment  = 1)
#'
#' fitSurv <- bdpsurvival(Surv(time, status) ~ historical + treatment,
#'                        data = example_surv_1arm)
#'
#' summary(fitSurv)
#'
#'
#' @rdname bdpsurvival
#' @import methods
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
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
           two_side      = TRUE){
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
           two_side      = TRUE){

  ### Check dataframe and ensure it has the correct column names
  namesData <- tolower(names(data))
  namesDiff <- setdiff(c("status", "time", "historical", "treatment"), namesData)
  if(length(namesDiff)>0){
    nDiff <- length(namesDiff)
    if(nDiff == 1){
      errorMsg <- paste0("Column ",
                         namesDiff,
                         " is missing from the input dataframe.")
      stop(errorMsg)
    } else if(nDiff>1){
      errorNames <- paste0(namesDiff, collapse = ", ")
      errorMsg <- paste0("Columns are missing from input dataframe: ",
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
    warning("Breaks vector includeded 0. The zero value was removed.")
  }



  ### Split the data on the breaks
  dataSplit <- survSplit(formula,
                         cut     = breaks,
                         start   = "start",
                         episode = "interval",
                         data    = data)


  ### If surv_time is null, replace with median time
  if(is.null(surv_time)){
    surv_time <- median(data$time)
  }


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


  ### Compute arm2, internal indicator of a two-arm trial
  if(nrow(S_c) == 0 & nrow(S0_c) == 0){
    arm2 <- FALSE
  } else{
    arm2 <- TRUE
  }


  ### Check inputs
  if(!arm2){
    if(nrow(S_t) == 0) stop("Current treatment data missing or input incorrectly.")
    if(nrow(S0_t) == 0) stop("Historical treatment data missing or input incorrectly.")
  } else if(arm2){
    if(nrow(S_t) == 0) stop("Current treatment data missing or input incorrectly.")
    if(nrow(S_c) == 0) stop("Current control data missing or input incorrectly.")
    if(nrow(S0_t) == 0 & nrow(S0_c)==0) stop("Historical data input incorrectly.")
  }

  posterior_treatment <- posterior_survival(
    S             = S_t,
    S0            = S0_t,
    alpha_max     = alpha_max[1],
    fix_alpha     = fix_alpha,
    a0            = a0,
    b0            = b0,
    surv_time     = surv_time,
    number_mcmc   = number_mcmc,
    weibull_shape = weibull_shape[1],
    weibull_scale = weibull_scale[1],
    two_side      = two_side,
    breaks        = breaks)

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
      breaks        = breaks)
  } else{
    posterior_control <- NULL
  }


  f1 <- final_survival(posterior_treatment = posterior_treatment,
                       posterior_control   = posterior_control,
                       arm2                = arm2,
                       surv_time           = surv_time,
                       breaks              = breaks)

  args1 <- list(S_t           = S_t,
                S_c           = S_c,
                S0_t          = S0_t,
                alpha_max     = alpha_max,
                fix_alpha     = fix_alpha,
                a0            = a0,
                b0            = b0,
                surv_time     = surv_time,
                number_mcmc   = number_mcmc,
                weibull_scale = weibull_scale,
                weibull_shape = weibull_shape,
                two_side      = two_side,
                breaks        = breaks)

  me <- list(posterior_treatment = posterior_treatment,
             posterior_control   = posterior_control,
             f1                  = f1,
             args1               = args1)

  class(me) <- "bdpsurvival"

  return(me)
})


################################################################################
### Helper functions
################################################################################
### Estimate discount function weight for prior data assuming survival outcome
### - Use approximation to the hazard ratio
discount_function_survival <- function(S, S0, alpha_max, fix_alpha, a0, b0,
                                       number_mcmc, weibull_shape, weibull_scale,
                                       two_side){

  ### Extract intervals and count number of intervals
  ### - Below, S_int should equal S0_int
  S_int  <- levels(S$interval)
  S0_int <- levels(S0$interval)
  nInt   <- length(S_int)
  n0Int  <- length(S0_int)

  interval <- NULL

  ### Compute posterior of hazard rate comparing current and historical
  a_post <- b_post <- numeric(nInt)
  a_post0 <- b_post0 <- numeric(n0Int)

  hazard_post_aug_t <- hazard_post_aug_t0 <- matrix(NA, number_mcmc, nInt)

  ### Compute posterior values
  for(i in 1:nInt){
    a_post[i] <- a0 + sum(subset(S, interval==S_int[i])$status)
    b_post[i] <- b0 + sum(subset(S, interval==S_int[i])$exposure)

    a_post0[i] <- a0 + sum(subset(S0, interval==S0_int[i])$status)
    b_post0[i] <- b0 + sum(subset(S0, interval==S0_int[i])$exposure)

    ### Add on a very small value to avoid underflow
    hazard_post_aug_t[,i]  <- rgamma(number_mcmc, a_post[i],  b_post[i])+1e-12
    hazard_post_aug_t0[,i] <- rgamma(number_mcmc, a_post0[i], b_post0[i])+1e-12
  }

  ### Weight historical data via (approximate) hazard ratio comparing
  ### current vs historical
  R0     <- log(hazard_post_aug_t0)-log(hazard_post_aug_t)
  V0     <- 1/apply(R0,2,var)
  logHR0 <- R0%*%V0/sum(V0)  #weighted average  of SE^2

  p_test <- mean(logHR0 > 0)   #larger is higher failure

  if(fix_alpha == TRUE){
    alpha_discount <- alpha_max
  } else{
    if (!two_side) {
      alpha_discount <- pweibull(p_test, shape=weibull_shape, scale=weibull_scale)*alpha_max
    } else if (two_side){
      p_test1    <- ifelse(p_test > 0.5, 1 - p_test, p_test)
      alpha_discount <- pweibull(p_test1, shape=weibull_shape, scale=weibull_scale)*alpha_max
    }
  }

  return(list(alpha_discount = alpha_discount,
              p_hat          = p_test,
              posterior_flat = hazard_post_aug_t,
              prior          = hazard_post_aug_t0))
}


### Posterior estimation for piecewise exponential distribution given
### alpha_loss.
### - For 1 arm (OPC), comparison is probability of survival at user input time
### - For 2 arm (RCT), comparison is hazard ratio of treatment vs control
posterior_augment_survival <- function(S, S0, alpha_discount, a0, b0,
                                       number_mcmc, surv_time, breaks){

  ### Extract intervals and count number of intervals
  ### - Below, S_int should equal S0_int
  S_int  <- levels(S$interval)
  S0_int <- levels(S0$interval)
  nInt   <- length(S_int)
  n0Int  <- length(S0_int)

  interval <- NULL

  ### Compute posterior of hazard rate comparing current and historical
  a_post <- b_post <- numeric(nInt)
  a_post0 <- b_post0 <- numeric(n0Int)

  hazard_post_aug_t <- matrix(NA, number_mcmc, nInt)

  for(i in 1:nInt){
    a_post0[i] <- a0 + sum(subset(S0, interval==S0_int[i])$status)
    b_post0[i] <- b0 + sum(subset(S0, interval==S0_int[i])$exposure)

    a_post[i] <- a0 + sum(subset(S, interval==S_int[i])$status)
    b_post[i] <- b0 + sum(subset(S, interval==S_int[i])$exposure)

    ### Add on a very small value to avoid underflow
    hazard_post_aug_t[,i]  <- rgamma(number_mcmc,
                                     a_post[i]+a_post0[i]*alpha_discount,
                                     b_post[i]+b_post0[i]*alpha_discount) + 1e-12
 }

  ### Posterior of survival time (statistic of interest for 1arm)
  pwe_cdf <- ppexp(q    = surv_time,
                   x    = hazard_post_aug_t,
                   cuts = c(0,breaks))

  surv_time_posterior <- 1-pwe_cdf

  return(list(cdf                 = pwe_cdf,
              surv_time_posterior = surv_time_posterior,
              posterior           = hazard_post_aug_t))
}


### Combine  loss function and posterior estimation into one function
posterior_survival <- function(S, S0, alpha_max, fix_alpha, a0, b0, surv_time,
                               number_mcmc, weibull_shape, weibull_scale,
                               two_side, breaks){

  alpha_discount <- discount_function_survival(S             = S,
                                               S0            = S0,
                                               alpha_max     = alpha_max,
                                               fix_alpha     = fix_alpha,
                                               a0            = a0,
                                               b0            = b0,
                                               number_mcmc   = number_mcmc,
                                               weibull_shape = weibull_shape,
                                               weibull_scale = weibull_scale,
                                               two_side      = two_side)

  posterior <- posterior_augment_survival(
    S              = S,
    S0             = S0,
    alpha_discount = alpha_discount$alpha_discount,
    a0             = a0,
    b0             = b0,
    number_mcmc    = number_mcmc,
    surv_time      = surv_time,
    breaks         = breaks)

    return(list(alpha_discount = alpha_discount$alpha_discount,
                p_hat          = alpha_discount$p_hat,
                posterior      = posterior,
                posterior_flat = alpha_discount$posterior_flat,
                prior          = alpha_discount$prior))
}

### Create final result class
final_survival <- function(posterior_treatment, posterior_control, arm2=FALSE,
                           surv_time, breaks){

  ### if(opc) plot survival prob.
  treatment_posterior      <- posterior_treatment$posterior$surv_time_posterior
  treatment_posterior_flat <- ppexp(surv_time, posterior_treatment$posterior_flat, c(0,breaks))
  treatment_prior          <- ppexp(surv_time, posterior_treatment$prior, c(0,breaks))

  density_post_treatment  <- density(treatment_posterior,
                                     adjust = 0.5)
  density_flat_treatment  <- density(treatment_posterior_flat,
                                     adjust = 0.5)
  density_prior_treatment <- density(treatment_prior,
                                     adjust = 0.5)

  if(!arm2){
    treatmentpost <- posterior_treatment$posterior$surv_time_posterior

    return(list(density_post_treatment  = density_post_treatment,
                density_flat_treatment  = density_flat_treatment,
                density_prior_treatment = density_prior_treatment,
                treatmentpost           = treatmentpost))
  } else if(!is.null(posterior_control) & arm2){
    ### Finalize below code
    density_post_control  <- density(posterior_control$posterior$posterior,
                                     adjust = 0.5)
    density_flat_control  <- density(posterior_control$posterior_flat,
                                     adjust = 0.5)
    density_prior_control <- density(posterior_control$prior,
                                     adjust = 0.5)

    ### Posterior hazard ratios at each interval
    R0     <- log(posterior_treatment$posterior$posterior)-log(posterior_control$posterior$posterior)
    V0     <- 1/apply(R0,2,var)
    logHR0 <- R0%*%V0/sum(V0)

    treatmentpost <- logHR0

    return(list(density_post_treatment  = density_post_treatment,
                density_flat_treatment  = density_flat_treatment,
                density_prior_treatment = density_prior_treatment,
                density_post_control    = density_post_control,
                density_flat_control    = density_flat_control,
                density_prior_control   = density_prior_control,
                treatmentpost           = treatmentpost))
  }
}


### Transform posterior hazard density list into a data frame
hazard_list_to_df <- function(hazard, starts, information_sources, group){
  if(!is.list(hazard)){
    return(hazard)
  }

  nS <- length(starts)
  options(warn=-1)
  df <- data.frame(information_sources = information_sources,
                   group = group,
                   start = starts[1],
                   x = hazard[[1]]$x,
                   y = hazard[[1]]$y)

  for(i in 2:nS){
    df1 <- data.frame(information_sources = information_sources,
                      group = group,
                      start = starts[i],
                      x = hazard[[i]]$x,
                      y = hazard[[i]]$y)
    df <- rbind(df, df1)
  }
  options(warn=0)

  return(df)
}
