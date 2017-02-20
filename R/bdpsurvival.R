#' bdpsurv
#'
#' bdpsurv
#'
#' @title bdpsurv: bdpsurv
#' @param formula an object of class "formula." Must have a survival object on the left side and exactly two inputs on the right side: treatment and historical. See ‘Details’ for more information.
#' @param data a data frame containing columns time, status, treatment, and historical.
#' @param breaks a vector of breaks used to compose the breaks of the piecewise exponential model.
#' @param a0 prior value for the gamma shape. Default is 1.
#' @param b0 prior value for the gamma rate. Default is 1.
#' @param surv_time survival time of interest for computing the the probability of survival for a type="1arm", i.e., an OPC trial. Default is median survival time.
#' @param type One of "1arm" or "2arm", denoting an OPC trial or a randomized control trial(RCT), respectively.
#' @param alpha_max Maximum weight the discount function can apply. Default is 1. For type="2arm", users may specify a vector of two values where the first value is used to weight the historical treatment group and the second value is used to weight the historical control group.
#' @param number_mcmc Number of Markov Chain Monte Carlo (MCMC) simulations. Default is 1e4.
#' @param weibull_shape Shape parameter of the Weibull discount function used to compute alpha, the weight parameter of the historical data. Default value is 3. For type="2arm", users may specify a vector of two values where the first value is used to estimate the weight of the historical treatment group and the second value is used to estimate the weight of the historical control group.
#' @param weibull_scale Scale parameter of the Weibull discount function used to compute alpha, the weight parameter of the historical data. Default value is 0.135. Two values have special treatment: 0 and Inf. For weibull_scale = 0, alpha is set to 0, i.e., no weight. For weibull_scale = Inf, alpha is set to 1, i.e., full weight. For type="2arm", users may specify a vector of two values where the first value is used to estimate the weight of the historical treatment group and the second value is used to estimate the weight of the historical control group.
#' @param two_side Indicator of two-sided test for the discount function. Default value is 1.
#'
#' @description
#'
#' @details
#'
#' @examples
#'
#' @rdname bdpsurvival
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
           type          = NULL,
           alpha_max     = 1,
           number_mcmc   = 10000,
           weibull_scale = 0.135,
           weibull_shape = 3,
           two_side      = 1){
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
           type          = NULL,
           alpha_max     = 1,
           number_mcmc   = 10000,
           weibull_scale = 0.135,
           weibull_shape = 3,
           two_side      = 1){

  if(type=="2arm"){
    return("Error: currently, only 1 arm (OPC) analyses are supported.")
  }


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

  ### If too many breaks, return error
  if(length(breaks) > 5){
    return("Error: currently only a maximum of 5 breaks are supported.")
  }


  ### Split the data on the breaks
  dataSplit <- survSplit(formula,
                         cut     = breaks,
                         start   = "start",
                         episode = "interval",
                         data    = data)


  ### If surv_time is null, replace with median time
  if(is.null(surv_time)){
    surv_time <- median(dataSplit$time)
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

  if(type=="1arm"){
    if(nrow(S_t) == 0) return("Error: current treatment data missing or input incorrectly.")
    if(nrow(S0_t) == 0) return("Error: historical treatment data missing or input incorrectly.")
  } else if(type=="2arm"){
    if(nrow(S_t) == 0) return("Error: current treatment data missing or input incorrectly.")
    if(nrow(S_c) == 0) return("Error: current control data missing or input incorrectly.")
    if(nrow(S0_t) == 0 & nrow(S0_c)) return("Error: historical data input incorrectly.")
  }

  posterior_treatment <- survival_posterior(
    S             = S_t,
    S0            = S0_t,
    alpha_max     = alpha_max[1],
    a0            = a0,
    b0            = b0,
    surv_time     = surv_time,
    number_mcmc   = number_mcmc,
    weibull_shape = weibull_shape[1],
    weibull_scale = weibull_scale[1],
    two_side      = two_side)

  if(type=="2arm"){
    posterior_control <- survival_posterior(
      S             = S_t,
      S0            = S0_t,
      alpha_max     = alpha_max[2],
      a0            = a0,
      b0            = b0,
      surv_time     = surv_time,
      number_mcmc   = number_mcmc,
      weibull_shape = weibull_shape[2],
      weibull_scale = weibull_scale[2],
      two_side      = two_side)
  }


  f1 <- final_survival(posterior_treatment = posterior_treatment)

  args1 <- list(S_t           = S_t,
                S_c           = S_c,
                S0_t          = S0_t,
                alpha_max     = alpha_max,
                a0            = a0,
                b0            = b0,
                surv_time     = surv_time,
                number_mcmc   = number_mcmc,
                weibull_scale = weibull_scale,
                weibull_shape = weibull_shape,
                two_side      = two_side)

  me <- list(posterior_treatment = posterior_treatment,
             f1                  = f1,
             args1               = args1)

  class(me) <- "bdpsurvival"

  return(me)
})


#' plot
#'
#' plot
#'
#' @title plot: plot
#' @param x bdpsurvival
#'
#' @examples
#'
#' @method plot
#'
#' @rdname plot
#' @export plot

setMethod("plot", signature(x = "bdpsurvival"), function(x){

  f                   <- x$f1
  posterior_treatment <- x$posterior_treatment
  two_side            <- x$args1$two_side
  starts              <- c(0,x$args1$breaks)

  D4 <- data.frame(information_sources = "Posterior",
                   group               = "Treatment",
                   y                   = f$density_post_treatment$y,
                   x                   = f$density_post_treatment$x)

  D5 <- data.frame(information_sources = "Current data",
                   group               = "Treatment",
                   y                   = f$density_flat_treatment$y,
                   x                   = f$density_flat_treatment$x)

  D6 <- data.frame(information_sources = "Prior",
                   group               = "Treatment",
                   y                   = f$density_prior_treatment$y,
                   x                   = f$density_prior_treatment$x)

  D <- as.data.frame(rbind(D4, D5, D6))

  D$information_sources <- factor(D$information_sources,
                                  levels = (c("Posterior", "Current data", "Prior")))

  D$start <- paste0("Interval start: ", D$start)

  post_typeplot <- ggplot(D, aes(x = x, y = y)) +
    geom_line(size = 2, aes(colour = information_sources, lty = information_sources)) +
    theme_bw() +
    facet_wrap(~group+start, ncol = 1, scale = "free") +
    ylab("Density (PDF)") +
    xlab("values")

  densityplot <- ggplot(subset(D, information_sources == "Posterior"), aes(x = x, y = y)) +
    geom_line(size = 2, aes(colour = group)) +
    ylab("Density (PDF)") +
    xlab("values") +
    theme_bw()


  if (two_side == 1) {
    p_value = seq(0, 1, , 100)
    p_value = ifelse(p_value > 0.5, 1 - p_value, p_value)
  }
  if (two_side == 0) {
    p_value = seq(0, 1, , 100)
  }

  Loss_function_treatment <- pweibull(p_value,
                                      shape = posterior_treatment$weibull_shape,
                                      scale = posterior_treatment$weibull_scale)

  D1 <- data.frame(group = "treatment", y = Loss_function_treatment, x = seq(0, 1, , 100))
  D2 <- data.frame(group = c("treatment"), pvalue = c(posterior_treatment$pvalue))

  lossfun_plot <- ggplot() +
    geom_line(data = D1, aes(y = y, x = x, colour = group), size = 1) +
    geom_vline(data = D2, aes(xintercept = pvalue, colour = group), lty = 2) +
    facet_wrap(~group, ncol = 1) +
    theme_bw() +
    ylab("Effective sample size for historical data") +
    xlab("Bayesian p-value (new vs historical data)")

  op <- par(ask=TRUE)
  plot(post_typeplot)
  plot(densityplot)
  plot(lossfun_plot)
  par(op)
})




#' summary
#'
#' summary
#'
#' @title summary: summary
#' @param object bdpsurvival
#'
#' @examples
#'
#' @method summary
#'
#' @rdname summary
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







################################################################################
### Helper functions
################################################################################
### Estimate discount function weight for prior data assuming survival outcome
### - Use approximation to the hazard ratio
discount_function_survival <- function(S, S0, alpha_max, a0, b0, number_mcmc,
                                       weibull_shape, weibull_scale, two_side){

  ### Extract intervals and count number of intervals
  ### - Below, S_int should equal S0_int
  S_int  <- levels(S$interval)
  S0_int <- levels(S0$interval)
  nInt   <- length(S_int)
  n0Int  <- length(S0_int)

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

  R0     <- log(hazard_post_aug_t0)-log(hazard_post_aug_t)
  V0     <- 1/apply(R0,2,var)
  logHR0 <- R0%*%V0/sum(V0)  #weighted average  of SE^2

  p_test <- mean(logHR0 > 0)   #larger is higher failure

  if(weibull_shape %in% c(0,Inf)){
    if(weibull_shape == 0){
      alpha_discount <- 0
    } else{
      alpha_discount <- 1
    }
  } else{
    if (two_side == 0) {
      alpha_discount <- pweibull(p_test, shape=weibull_shape, scale=weibull_scale)*alpha_max
    } else if (two_side == 1){
      p_test1    <- ifelse(p_test > 0.5, 1 - p_test, p_test)
      alpha_discount <- pweibull(p_test1, shape=weibull_shape, scale=weibull_scale)*alpha_max
    }
  }

  return(list(alpha_discount = alpha_discount,
              pvalue         = p_test,
              posterior_flat = hazard_post_aug_t,
              prior          = hazard_post_aug_t0))
}


### Posterior estimation for piecewise exponential distribution given
### alpha_loss. Comparison is probability of survival at user input time
posterior_augment_survival <- function(S, S0, alpha_discount, a0, b0,
                                       number_mcmc, surv_time){

  ### Extract intervals and count number of intervals
  ### - Below, S_int should equal S0_int
  S_int  <- levels(S$interval)
  S0_int <- levels(S0$interval)
  nInt   <- length(S_int)
  n0Int  <- length(S0_int)

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

  pwe_cdf <- ppexpV(surv_time = surv_time,
                    hazard    = hazard_post_aug_t,
                    breaks    = breaks)

  surv_time_posterior <- 1-pwe_cdf

  return(list(cdf                 = pwe_cdf,
              surv_time_posterior = surv_time_posterior,
              posterior           = hazard_post_aug_t))
}


### Combine  loss function and posterior estimation into one function
survival_posterior <- function(S, S0, alpha_max, a0, b0, surv_time,
                               number_mcmc, weibull_shape, weibull_scale,
                               two_side){

  alpha_discount <- discount_function_survival(S             = S,
                                               S0            = S0,
                                               alpha_max     = alpha_max,
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
    surv_time      = surv_time)

    return(list(alpha_discount = alpha_discount$alpha_discount,
                pvalue         = alpha_discount$pvalue,
                posterior      = posterior,
                posterior_flat = alpha_discount$posterior_flat,
                prior          = alpha_discount$prior,
                weibull_scale  = weibull_scale,
                weibull_shape  = weibull_shape,
                S              = S,
                S0             = S0))
}

### Create final result class
final_survival <- function(posterior_treatment){
  density_post_treatment  <- density(posterior_treatment$Survival_posterior,
                                   adjust = 0.5)
  density_flat_treatment  <- density(posterior_treatment$posterior_flat,
                                   adjust = 0.5)
  density_prior_treatment <- density(posterior_treatment$prior,
                                   adjust = 0.5)

  teatmentpost <- posterior_treatment$posterior$surv_time_posterior

  return(list(density_post_treatment  = density_post_treatment,
              density_flat_treatment  = density_flat_treatment,
              density_prior_treatment = density_prior_treatment,
              teatmentpost            = teatmentpost))
}



### Piecewise exponential cdf
ppexp <- function (q, rate = 1, t = 0){
  q[q < 0] <- 0
  ind <- rowSums(outer(q, t, ">="))
  ret <- pexp(q - t[ind], rate[ind])
  mi <- min(length(t), max(ind))
  if (length(t) > 1) {
    dt <- t[-1] - t[-mi]
    pe <- pexp(dt, rate[-mi])
    cp <- c(1, cumprod(1 - pe))
    ret <- c(0, cumsum(cp[-length(cp)] * pe))[ind] + ret * cp[ind]
  }
  ret
}

### Vectorized ppexp function over rates
ppexpV <- function(q, rate, t){
  nR  <- nrow(rate)
  ret <- apply(rate, 1, ppexp, q=q, t=t)
  ret
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
