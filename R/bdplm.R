#' @title Bayesian Discount Prior: Two-Arm Linear Regression
#' @description \code{bdplm} is used to estimate the treatment effect
#'   in the presence of covariates using the regression analysis
#'   implementation of the Bayesian discount prior. This function is a
#'   barebones implementation and does not support \code{summary},
#'   \code{plot}, and \code{print} methods.
#' @param formula an object of class "formula." See "Details" for
#'   more information, including specfication of historical and treatment
#'   data indicators.
#' @param data an optional data frame, list or environment
#'   (or object coercible by as.data.frame to a data frame)
#'   containing the variables in the model. If not found in data,
#'   the variables are taken from environment(formula), typically
#'   the environment from which bdplm is called.
#' @param prior_mean historical mean for the coefficients, specified as a
#'   a vector. The length of \code{prior_mean} must be equal to the number
#'   of adjusting covariates + 2 (ignore the intercept in this case). See
#'   "Details" for more information.
#' @param prior_sigma prior standard deviation for the coefficient, specified
#'   as a vector. The length of \code{prior_sigma} must be equal to the number
#'   of adjusting covariates + 2 (ignore the intercept in this case). See
#'   "Details" for more information.
#' @param alpha_max scalar. Maximum weight the discount function can apply.
#'   Default is 1. Users may specify a vector of two values where the first
#'   value is used to weight the historical treatment group and
#'   the second value is used to weight the historical control group.
#' @param fix_alpha logical. Fix alpha at alpha_max? Default value is FALSE.
#' @param number_mcmc_alpha scalar. Number of Monte Carlo
#'   simulations for estimating the historical data weight. Default is 10000.
#' @param number_mcmc_sigmagrid scalar. Grid size for computing sigma.
#'   Default is 1000. See "Details" for more information.
#' @param number_mcmc_sigma scalar. Number of Monte Carlo simulations for
#'   estimating sigma. Default is 1000. See "Details" for more information.
#' @param number_mcmc_beta scalar. Number of Monte Carlo simulations for
#'   estimating beta, the vector of regression coefficients. Default is 10000.
#' @param weibull_shape scalar. Shape parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 3. Users may specify a vector of two values
#'   where the first value is used to estimate the weight of the historical
#'   treatment group and the second value is used to estimate the weight of the
#'   historical control group.
#' @param weibull_scale scalar. Scale parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 0.135. Users may specify a vector of two
#'   values where the first value is used to estimate the weight of the
#'   historical treatment group and the second value is used to estimate the
#'   weight of the historical control group.
#' @param method character. Analysis method with respect to estimation of the
#'   weight paramter alpha. Default value "\code{fixed}" estimates alpha once
#'   and holds it fixed throughout the analysis. Alternative method
#'   "\code{mc}" estimates alpha for each Monte Carlo iteration. Currently, only
#'   the default method "\code{fixed}" is supported.
#' @details \code{bdplm} uses a two-stage approach for determining the
#'   strength of historical data in estimation of an adjusted mean or covariate
#'   effect. In the first stage, a Weibull distribution function is used as a
#'   \emph{discount function} that defines the maximum strength of the
#'   historical data (via \code{weibull_shape}, \code{weibull_scale}, and
#'   \code{alpha_max}) and discounts based on disagreement with the current data.
#'   Disagreement between current and historical data is determined by stochastically
#'   comparing the respective posterior distributions under noninformative priors.
#'   Here with a two-arm regression analysis, the comparison is the
#'   proability (\code{p}) that the covariate effect of an historical data indicator is
#'   significantly different from zero. The comparison metric \code{p} is then
#'   input into the Weibull discount function and the final strength of the
#'   historical data is returned (alpha).
#'
#'   In the second stage, posterior estimation is performed where the discount
#'   function parameter, \code{alpha}, is used as a fixed value for all posterior
#'   estimation procedures.
#'
#'   At minimum, the formula/data must include an intercept and each of historical
#'   and treatment columns. Any covariates can be included as well. Example usage
#'   that includes the intercept by default could be:
#'   \code{y ~ hisorical+treatment+baseline}.
#'
#'   In this implementation, current and historical data must be present for both
#'   treatment and control groups. The covariate prior parameters \code{prior_mean}
#'   and \code{prior_sigma} must be input as a vector with elements in the order
#'   given by the formula.
#'
#' @return \code{bdplm} returns an object of class "bdplm".
#'
#' An object of class "\code{bdplm}" is a list containing at least
#' the following components:
#' \describe{
#'    \itemize{
#'      \item{\code{posterior}}{
#'        data frame. The posterior draws of the covariates, the intercept, and
#'        the treatment effect. The grid of sigma values are included.}
#'      \item{\code{alpha_discount}}{
#'        numeric. The posterior probability of the stochastic comparison
#'        between the current and historical data.}
#'    }
#' }
#'
#' @examples
#' # Simulate  data
#' n_t  <- 100
#' n_c  <- 100
#' n_t0 <- 250
#' n_c0 <- 250
#' x          <- rnorm(n_t+n_c+n_t0+n_c0, 34, 5)
#' treatment  <- c(rep(1, n_t+n_t0), rep(0, n_c+n_c0))
#' historical <- c(rep(0, n_t), rep(1,n_t0), rep(0, n_c), rep(1,n_c0))
#'
#' Y <- treatment + 1000*historical + x*3.5 + rnorm(n_t+n_c+n_t0+n_c0)
#' df <- data.frame(Y=Y, treatment=treatment, historical=historical, x=x)
#'
#' fit <- bdplm(Y ~ treatment+historical+x, data=df,
#'              prior_mean  = rep(0,3),
#'              prior_sigma = rep(1e4,3))
#'
#'
#' @rdname bdplm
#' @import methods
#' @importFrom stats density is.empty.model median model.offset model.response pweibull pnorm quantile rbeta rgamma rnorm var vcov contrasts<- dt gaussian lm.fit model.frame model.matrix.default offset rchisq terms terms.formula coefficients pchisq
#' @aliases bdplm,ANY-method
#' @useDynLib bayesDP
#' @export bdplm
bdplm <- setClass("bdplm", slots = c(posterior_treatment = "list",
                                     posterior_control = "list",
                                     args1 = "list"))
setGeneric("bdplm",
  function(formula                   = formula,
           data                      = data,
           prior_mean                = 0,
           prior_sigma               = 1,
           number_mcmc_alpha         = 10000,
           number_mcmc_sigmagrid     = 1000,
           number_mcmc_sigma         = 100,
           number_mcmc_beta          = 10000,
           alpha_max                 = 1,
           fix_alpha                 = FALSE,
           weibull_scale             = 0.135,
           weibull_shape             = 3,
           method                    = "fixed"){
             standardGeneric("bdplm")
           })

setMethod("bdplm",
  signature(),
  function(formula                   = formula,
           data                      = data,
           prior_mean                = 0,
           prior_sigma               = 1,
           number_mcmc_alpha         = 10000,
           number_mcmc_sigmagrid     = 1000,
           number_mcmc_sigma         = 100,
           number_mcmc_beta          = 10000,
           alpha_max                 = 1,
           fix_alpha                 = FALSE,
           weibull_scale             = 0.135,
           weibull_shape             = 3,
           method                    = "fixed"){

  ### Check validity of family input
  call <- match.call()
  ### Pull data from the environment if data input is missing
  if (missing(data)) {
    data <- environment(formula)
  }

  ### Check method
  if(method != "fixed"){
    stop("Only method = 'fixed' is currently supported.")
  }

  ### Place data into X and Y objects
  mf <- match.call(expand.dots = FALSE)
  m  <- match(c("formula", "data", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$na.action <- NULL
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) {
      names(Y) <- nm
    }
  }

  X <- if (!is.empty.model(mt)) {
    model.matrixBayes(object = mt, data = data, contrasts.arg = NULL,
                      keep.order = TRUE, drop.baseline = TRUE)
  } else {
    matrix(, NROW(Y), 0L)
  }

  ### Alter intercept column name, if present
  imatch <- match("(Intercept)", colnames(X))
  if(!is.na(imatch)) colnames(X)[imatch] <- "intercept"

  ### Assign one-arm or two-arm analysis. Check validity of inputs.
  # Create indicator of whether each column is present
  cmatch <- match(c("intercept", "historical", "treatment") , colnames(X))
  cmatch <- !is.na(cmatch)  ### == TRUE == present

  if(!all(cmatch)){
    stop("Data is input incorrectly. Intercept, historical, and treatment must all be present.")
  }


  ### Count levels of historical data and ensure 0 and 1 are present
  hist_levels <- levels(as.factor(X[,"historical"]))

  # Check that the historical levels are 0 and 1
  if(!(all(hist_levels %in% c(0,1)))){
    stop("Historical input has wrong levels. Values should be 0 and 1.")
  }


  ### Count levels of treatment data and sure 0 and 1 are present
  trt_levels <- levels(as.factor(X[,"treatment"]))

  # Check that the levels are 0 and/or 1
  if(!(all(trt_levels %in% c(0,1)))){
    stop("Treatment input has wrong levels. Values should be 0 and 1.")
  }


  historical <- NULL
  treatment  <- NULL
  intercept  <- NULL

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
  ### Create a master dataframe
  df <- data.frame(Y=Y, data.frame(X))

  ### Split data into separate treatment & control dataframes
  df_t <- subset(df, treatment == 1)
  df_c <- subset(df, treatment == 0)

  ### Also create current data dataframe
  df_current <- subset(df, historical == 0, select=-intercept)

  ### Count number of covariates
  names_df <- names(df)
  covs_df  <- names_df[!(names_df %in% c("Y", "treatment", "historical", "intercept"))]
  n_covs   <- length(covs_df)


  ##############################################################################
  # Estimate discount weights for each of the treatment and control arms
  ##############################################################################
  discount_treatment <- discount_lm(df                = df_t,
                                    alpha_max         = alpha_max[1],
                                    fix_alpha         = fix_alpha,
                                    number_mcmc_alpha = number_mcmc_alpha,
                                    weibull_shape     = weibull_shape[1],
                                    weibull_scale     = weibull_scale[1])

  discount_control <- discount_lm(df                = df_c,
                                  alpha_max         = alpha_max[2],
                                  fix_alpha         = fix_alpha,
                                  number_mcmc_alpha = number_mcmc_alpha,
                                  weibull_shape     = weibull_shape[2],
                                  weibull_scale     = weibull_scale[2])


  ##############################################################################
  # Estimate augmented treatment effect
  ##############################################################################
  ### Compute prior terms
  tau2   <- prior_sigma^2
  mu0    <- prior_mean

  ### Extract alpha0, append "zero" weight for the covariate effect(s)
  alpha0 <- c(discount_treatment$alpha_discount + 1e-12,
              discount_control$alpha_discount + 1e-12,
              rep(1e-12, n_covs))

  ### Calculate constants from current data
  df_current$control <- 1-df_current$treatment

  X     <- df_current[,c("treatment", "control",covs_df)]
  X     <- as.matrix(X)
  y     <- df_current$Y
  ystar <- c(y, mu0)
  Xstar <- rbind(X, diag(length(mu0)))
  XtX   <- t(X)%*%X
  Xty   <- t(X)%*%y
  SigmaBetaInv <- diag(alpha0/tau2)

  ### Grid search of sigma2
  ### - Find grid limits via wls
  W     <- c(rep(1,length(y)), alpha0/tau2)
  lmfit <- lm(ystar~Xstar-1, weights=W)
  summ  <- summary(lmfit)
  a     <- summ$df[2]
  b     <- summ$sigma^2
  lower <- 1/qgamma(.9999999999,a/2,(a*b)/2)
  upper <- 1/qgamma(.0000000001,a/2,(a*b)/2)

  grid_sigma2 <- seq(lower,upper,length.out=number_mcmc_sigmagrid)

  ### Sample candidate values of sigma2
  sigma2candidates <- sigma2marginal(n           = number_mcmc_sigmagrid,
                                    grid         = grid_sigma2,
                                    XtX          = XtX,
                                    SigmaBetaInv = SigmaBetaInv,
                                    Xstar        = Xstar,
                                    Xty          = Xty,
                                    mu0          = mu0,
                                    ystar        = ystar)

  ### Normalize the marginal posteriors (log-likelihoods) and exponentiate
  logL  <- sigma2candidates$logL
  logL  <- logL[is.finite(logL)]
  normL <- logL[which.min(abs(logL))]
  L     <- exp(logL - normL)

  ### Sample with replacement from marginal posterior density of sigma2
  sigma2_accept <- sample(x       = sigma2candidates$sigma2,
                          size    = number_mcmc_sigma,
                          replace = TRUE,
                          prob    = L)


  ### Draw samples of the covariate vector beta
  ### n_beta_samples is ceiling(number_mcmc_beta/number_mcmc_sigma2)
  n_beta_samples <- ceiling(number_mcmc_beta/number_mcmc_sigma)

  beta_samples <- betaRegSampler(sigma2_accept, XtX, SigmaBetaInv, mu0,
                                 Xty, n_beta_samples)

  beta_samples <- data.frame(beta_samples)
  names(beta_samples) <- c("treatment", "control", covs_df, "sigma")

  ### Estimate posterior of intercept and treatment effect
  beta_samples$intercept <- beta_samples$control
  beta_samples$treatment <- beta_samples$treatment-beta_samples$control
  beta_samples$sigma     <- sqrt(beta_samples$sigma)

  ### Format alpha_discount values
  alpha_discount <- data.frame(treatment = discount_treatment$alpha_discount,
                               control   = discount_control$alpha_discount)

  ### Format estimates
  estimates <- list()
  estimates$coefficients <- data.frame(t(colMeans(posterior)))
  estimates$se           <- data.frame(t(apply(posterior,2,sd)))


  me <- list(posterior      = beta_samples,
             alpha_discount = alpha_discount,
             estimates      = estimates)

  class(me) <- "bdplm"
  return(me)
})




################################################################################
# Linear Regression discount weight estimation
# 1) Estimate the discount function
#    - Test that the historical effect is different from zero
#    - Estimate alpha based on the above comparison
################################################################################
discount_lm <- function(df, alpha_max, fix_alpha,
                         number_mcmc_alpha,
                         weibull_shape, weibull_scale){

  # Create formula
  cnames <- names(df)
  cnames <- cnames[!(cnames %in% c("Y", "intercept", "historical","treatment"))]
  cnames <- c("historical", cnames)
  f      <- paste0("Y~",paste0(cnames,collapse="+"))

  ### Get "flat" fit
  lm_fit   <- lm(f, data=df)
  lm_summ  <- summary(lm_fit)

  ### Monte Carlo simulations of error variance
  a        <- lm_summ$df[2]
  b        <- lm_summ$sigma^2
  sigma2   <- 1/rgamma(number_mcmc_alpha,a/2,(a*b)/2)

  ### Monte Carlo simulations of historical effect
  sigma2_beta <- lm_summ$cov.unscaled[2,2]*sigma2
  beta        <- rnorm(number_mcmc_alpha,lm_summ$coefficients[2,1],sqrt(sigma2_beta))

  p_hat  <- mean(beta>0)
  p_hat  <- ifelse(p_hat > 0.5, 1 - p_hat, p_hat)

  if(fix_alpha){
    alpha0 <- alpha_max
  } else{
    alpha0 <- pweibull(p_hat, shape = weibull_shape, scale = weibull_scale)*alpha_max
  }

  res <- list(p_hat                     = p_hat,
              alpha_discount            = alpha0)
  res
}









model.matrixBayes <- function(object, data=environment(object), contrasts.arg=NULL,
                              xlev=NULL, keep.order=FALSE, drop.baseline=FALSE,...){

  t <- if( missing( data ) ) {
    terms( object )
  } else{
    terms.formula(object, data = data, keep.order=keep.order)
  }

  attr(t, "intercept") <- attr(object, "intercept")
  if (is.null(attr(data, "terms"))){
    data <- model.frame(object, data, xlev=xlev)
  } else {
    reorder <- match(sapply(attr(t,"variables"), deparse, width.cutoff=500)[-1], names(data))
    if (anyNA(reorder)) {
      stop( "model frame and formula mismatch in model.matrix()" )
    }
    if(!identical(reorder, seq_len(ncol(data)))) {
      data <- data[,reorder, drop = FALSE]
    }
  }

  int <- attr(t, "response")
  if(length(data)) {      # otherwise no rhs terms, so skip all this
    if (drop.baseline){
      contr.funs <- as.character(getOption("contrasts"))
    }else{
      contr.funs <- as.character(list("contr.bayes.unordered", "contr.bayes.ordered"))
    }

    namD <- names(data)

    ## turn  character columns into factors
    for(i in namD)
      if(is.character( data[[i]] ) ) {
        data[[i]] <- factor(data[[i]])
        warning( gettextf( "variable '%s' converted to a factor", i ), domain = NA)
      }
    isF <- vapply(data, function(x) is.factor(x) || is.logical(x), NA)
    isF[int] <- FALSE
    isOF <- vapply(data, is.ordered, NA)
    for( nn in namD[isF] )            # drop response
      if( is.null( attr( data[[nn]], "contrasts" ) ) ) {
        contrasts( data[[nn]] ) <- contr.funs[1 + isOF[nn]]
      }

    ## it might be safer to have numerical contrasts:
    ##    get(contr.funs[1 + isOF[nn]])(nlevels(data[[nn]]))
    if ( !is.null( contrasts.arg ) && is.list( contrasts.arg ) ) {
      if ( is.null( namC <- names( contrasts.arg ) ) ) {
        stop( "invalid 'contrasts.arg' argument" )
      }
      for (nn in namC) {
        if ( is.na( ni <- match( nn, namD ) ) ) {
          warning( gettextf( "variable '%s' is absent, its contrast will be ignored", nn ), domain = NA )
        }
        else {
          ca <- contrasts.arg[[nn]]
          if( is.matrix( ca ) ) {
            contrasts( data[[ni]], ncol( ca ) ) <- ca
          }
          else {
            contrasts( data[[ni]] ) <- contrasts.arg[[nn]]
          }
        }
      }
    }
  } else {
    isF  <-  FALSE
    data <- data.frame(x=rep(0, nrow(data)))
  }
  ans  <- model.matrix.default(object=t, data=data)
  cons <- if(any(isF)){
    lapply( data[isF], function(x) attr( x,  "contrasts") )
  }else { NULL }
  attr(ans, "contrasts" ) <- cons
  ans
}
