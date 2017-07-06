#' @title Bayesian Discount Prior: Regression Analysis
#' @description \code{bdpregression} is used to estimate the adjusted intercept
#'   (single arm trial; OPC) or treatment effect (two-arm trial; RCT) for
#'   data using the regression analysis implementation of the
#'   Bayesian discount prior.
#' @param formula an object of class "formula." See "Details" for
#'   more information, including specfication of historical and treatment
#'   data status.
#' @param family 	a description of the error distribution and link
#'   function to be used in the model. For bdpregression this can
#'   be a character string naming a family function, a family
#'   function or the result of a call to a family function.
#' @param data an optional data frame, list or environment
#'   (or object coercible by as.data.frame to a data frame)
#'   containing the variables in the model. If not found in data,
#'   the variables are taken from environment(formula), typically
#'   the environment from which bdpregression is called.
#' @param prior_mean prior mean for the coefficients; default is 0.
#'   Can be a vector of length equal to the number of covariates
#'   (not counting the intercept, if any). If it is a scalar, it is
#'   expanded to the length of this vector.
#' @param prior_scale prior scale for the coefficients: default is NULL;
#'   if NULL, for a logistic regression model, prior_scale is 2.5; for a
#'    probit model, prior scale is 2.5*1.6. Can be a vector of length equal
#'    to the number of predictors (not counting the intercept, if any). If
#'    it is a scalar, it is expanded to the length of this vector.
#' @param prior_df prior degrees of freedom for the coefficients. For
#'    t distribution default is 1 (Cauchy). Set to Inf to get normal prior
#'    distributions. Can be a vector of length equal to the number of
#'    predictors (not counting the intercept, if any). If it is a scalar,
#'    it is expanded to the length of this vector.
#' @param prior_mean_for_intercept prior mean for the intercept: default
#'    is 0.
#' @param prior_scale_for_intercept prior scale for the intercept:
#'    default is NULL; for a logit model, prior scale for intercept is 10;
#'    for probit model, prior scale for intercept is rescaled as 10*1.6.
#' @param prior_df_for_intercept prior degrees of freedom for the
#'    intercept: default is 1.
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
#' @param method character. Analysis method with respect to estimation of the weight
#'   paramter alpha. Default value "\code{fixed}" estimates alpha once and holds it fixed
#'   throughout the analysis. Alternative method "\code{mc}" estimates alpha for each
#'   Monte Carlo iteration. Currently, only the default method "\code{fixed}" is
#'   supported.
#' @details \code{bdpregression} uses a two-stage approach for determining the
#'   strength of historical data in estimation of an adjusted mean or covariate effect.
#'   In the first stage, a Weibull distribution function is used as a
#'   \emph{discount function} that defines the maximum strength of the
#'   historical data (via \code{weibull_shape}, \code{weibull_scale}, and
#'   \code{alpha_max}) and discounts based on disagreement with the current data.
#'   Disagreement between current and historical data is determined by stochastically
#'   comparing the respective posterior distributions under noninformative priors.
#'   With a single arm regression analysis, the comparison is the
#'   proability (\code{p}) that the current adjusted intercept is less than the
#'   historical adjusted intercept. The comparison metric \code{p} is then
#'   input into the Weibull discount function and the final strength of the
#'   historical data is returned (alpha).
#'
#'   In the second stage, posterior estimation is performed where the discount
#'   function parameter, \code{alpha}, is used as a fixed value for all posterior
#'   estimation procedures.
#'
#'  Two-arm (RCT) analyses are not currently available with this release.
#'
#' @return \code{bdpregression} returns an object of class "bdpregression".
#'   The functions \code{summary} and \code{print} are used to obtain and
#'   print a summary of the results, including user inputs. The \code{plot}
#'   function displays visual outputs as well.
#'
#' An object of class "\code{bdpregression}" is a list containing at least
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
#'      \item{\code{posterior_regression}}{
#'        list. Contains results for the augmented regression analysis. Entries
#'        are similar to the output of \code{glm} and \code{bayesglm} (from the
#'         \code{arm} package).}
#'      \item{\code{posterior_flat_regression}}{
#'        list. Contains entries similar to \code{posterior_regression} corresponding
#'        to estimates of the unweighted current data.}
#'      \item{\code{prior_regression}}{
#`        list. Contains entries similar to \code{posterior_regression} corresponding
#'        to estimates of the historical data.}
#'    }
#'  \item{\code{args1}}{
#'    list. Entries contain user inputs. In addition, the following elements
#'    are ouput:}
#'    \itemize{
#'      \item{\code{df_t} and \code{df_c}}{
#'        dataframe. Input data parsed into internally used treatment and control
#'        data frames.
#'      }
#'      \item{\code{arm2}}{
#'        logical. Used internally to indicate one-arm or two-arm analysis.
#'      }
#'   }
#' }
#'
#' @examples
#' # One-arm trial (OPC) example - linear regression
#' # Simulate regression data for a single arm (OPC) trial
#' set.seed(3408)
#' historical <- c(rep(1,50),rep(0,50))
#' x1         <- c(rnorm(50), rnorm(50))
#' y          <- c(rnorm(50), rnorm(50)+0.2)
#'
#' fit1 <- bdpregression(y ~ x1 + historical)
#' summary(fit1)
#' print(fit1)
#' plot(fit1, type="discount")
#'
#'
#' # One-arm trial (OPC) example - logistic regression
#' set.seed(3408)
#' historical <- c(rep(1,100),rep(0,100))
#' x1         <- c(rep(0,50), rep(1,50), rep(0,50), rep(1,50))
#' y          <- rbinom(200,1,plogis(1 + 0.5*x1 + 0.1*historical))
#'
#' fit2 <- bdpregression(y ~ x1 + historical,
#'                       family = binomial)
#' print(fit2)
#'
#'
#' # One-arm trial (OPC) example - Poisson regression
#' set.seed(3408)
#' historical <- c(rep(1,100),rep(0,100))
#' x1         <- c(rep(0,50), rep(1,50), rep(0,50), rep(1,50))
#' y          <- rpois(200,exp(1 + 0.5*x1 + 0.5*historical))
#'
#' fit3 <- bdpregression(y ~ x1 + historical,
#'                       family = poisson)
#' summary(fit3)
#'
#' # Refit the Poisson data ignoring historical
#' fit4 <- bdpregression(y ~ x1, family = poisson)
#' summary(fit4)
#'
#'
#' # Place data in a dataframe and carry out linear regression
#' set.seed(3408)
#' df <- data.frame(historical = c(rep(1,50),rep(0,50)),
#'                  x1         = c(rnorm(50), rnorm(50)),
#'                  y          = c(rnorm(50), rnorm(50)+0.2))
#'
#' fit5 <- bdpregression(y ~ x1 + historical, data=df)
#' summary(fit5)
#'
#'
#' # Two-arm trials are not yet implemented.
#'
#'
#' @rdname bdpregression
#' @import methods
#' @importFrom stats density is.empty.model median model.offset model.response pweibull pnorm quantile rbeta rgamma rnorm var vcov contrasts<- dt gaussian lm.fit model.frame model.matrix.default offset rchisq terms terms.formula coefficients pchisq
#' @aliases bdpregression,ANY-method
#' @export bdpregression
bdpregression <- setClass("bdpregression", slots = c(posterior_treatment = "list",
                                                     posterior_control = "list",
                                                     args1 = "list"))
setGeneric("bdpregression",
  function(formula                   = formula,
           family                    = "gaussian",
           data                      = data,
           prior_mean                = 0,
           prior_scale               = NULL,
           prior_df                  = 1,
           prior_mean_for_intercept  = 0,
           prior_scale_for_intercept = NULL,
           prior_df_for_intercept    = 1,
           alpha_max                 = 1,
           fix_alpha                 = FALSE,
           number_mcmc               = 10000,
           weibull_scale             = 0.135,
           weibull_shape             = 3,
           two_side                  = TRUE,
           method                    = "fixed"){
             standardGeneric("bdpregression")
           })

setMethod("bdpregression",
  signature(),
  function(formula                   = formula,
           family                    = "gaussian",
           data                      = data,
           prior_mean                = 0,
           prior_scale               = NULL,
           prior_df                  = 1,
           prior_mean_for_intercept  = 0,
           prior_scale_for_intercept = NULL,
           prior_df_for_intercept    = 1,
           alpha_max                 = 1,
           fix_alpha                 = FALSE,
           number_mcmc               = 10000,
           weibull_scale             = 0.135,
           weibull_shape             = 3,
           two_side                  = TRUE,
           method                    = "fixed"){

  ### Check validity of family input
  call <- match.call()
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

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

  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of\n            observations)",
                    length(offset), NROW(Y)), domain = NA)
  } else{
    offset <- numeric(NROW(Y))
  }


  ### Alter intercept column name, if present
  imatch <- match("(Intercept)", colnames(X))
  if(!is.na(imatch)) colnames(X)[imatch] <- "intercept"


  ### Assign one-arm or two-arm analysis. Check validity of inputs.
  # Create indicator of whether each column is present
  cmatch <- match(c("intercept", "historical", "treatment") , colnames(X))
  cmatch <- !is.na(cmatch)  ### == TRUE == present

  ### If historical present, count levels
  if(cmatch[2]){
    hist_levels <- levels(as.factor(X[,"historical"]))

    # Check that the levels are 0 and/or 1
    if(any(!(hist_levels %in% c(0,1)))){
      stop("Historical input has wrong levels. Values should be 0 and/or 1.")
    }

    n_hist_levels <- length(hist_levels)
  }

  ### If treatment present, count levels
  if(cmatch[3]){
    trt_levels <- levels(as.factor(X[,"treatment"]))

    # Check that the levels are 0 and/or 1
    if(any(!(trt_levels %in% c(0,1)))){
      stop("Treatment input has wrong levels. Values should be 0 and/or 1.")
    }

    n_trt_levels <- length(trt_levels)
  }



  if(!cmatch[1] & cmatch[2] & !cmatch[3]){
    ### No intercept and no treatment indicator: cannot analyze as
    ### one-arm (or two-arm)
    stop("Data input incorrectly. Intercept and (optional) treatment
          covariate are missing. Cannot determine if one-arm or two-arm
          analysis.")
  } else if(all(!cmatch)){
    ### No intercept, historical, or treatment: cannot figure out what kind of analysis
    ### (one arm makes no sense here since intercept must be present)
    stop("Data input incorrectly. Intercept, historical covariate, and (optional) treatment
          covariate are missing. Cannot determine if one-arm or two-arm analysis.")
  } else if(all(cmatch)){
    if(n_trt_levels == 2){
      arm2 <- TRUE
    } else if(n_trt_levels == 1){
      arm2 <- FALSE
    }
  } else if(!cmatch[1] & cmatch[2] & cmatch[3]){
    if(n_trt_levels == 2){
      arm2 <- TRUE
    } else if(n_trt_levels == 1){
      arm2 <- FALSE
    }
  } else if(!cmatch[1] & !cmatch[2] & cmatch[3]){
    if(n_trt_levels == 2){
      arm2 <- TRUE
    } else if(n_trt_levels == 1){
      arm2 <- FALSE
    }
  } else if(cmatch[1] & !cmatch[2] & cmatch[3]){
    if(n_trt_levels == 2){
      arm2 <- TRUE
    } else if(n_trt_levels == 1){
      arm2 <- FALSE
    }
  } else if(cmatch[1] & !cmatch[2] & !cmatch[3]){
    arm2 <- FALSE
  } else if(cmatch[1] & cmatch[2] & !cmatch[3]){
    arm2 <- FALSE
  }

  historical <- NULL
  treatment  <- NULL

  if(arm2) stop("Two arm trials are not currently supported.")


  ### Check that (if present) historical and (if present) treatment columns
  ### have exactly 1-2 levels
  if(cmatch[2]){
    if(!(length(levels(as.factor(X[,"historical"]))) %in% c(1,2))){
      stop("Historical column does not have 1 or 2 unique values.
            Convert to binary indicator or fix inputs.")
    }
  }

  if(arm2){
    if(!(length(levels(as.factor(X[,"treatment"]))) %in% c(1,2))){
      stop("Treatment column does not have 1 or 2 unique values.
           Convert to binary indicator or fix inputs.")
    }
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
  ### Create a master dataframe
  df <- data.frame(Y=Y, data.frame(X), offset=offset)

  ### Split data into separate treatment & control dataframes
  if(!arm2){
    df_t <- df
    df_c <- NULL
  } else if(arm2){
    df_t <- subset(df, treatment == 1)
    df_c <- subset(df, treatment == 0)
  }


  ##############################################################################
  # Posterior analysis
  ##############################################################################
  posterior_treatment <- posterior_regression(
    df                        = df_t,
    family                    = family,
    alpha_max                 = alpha_max[1],
    fix_alpha                 = fix_alpha,
    prior_mean                = prior_mean,
    prior_scale               = prior_scale,
    prior_df                  = prior_df,
    prior_mean_for_intercept  = prior_mean_for_intercept,
    prior_scale_for_intercept = prior_scale_for_intercept,
    prior_df_for_intercept    = prior_df_for_intercept,
    number_mcmc               = number_mcmc,
    weibull_shape             = weibull_shape[1],
    weibull_scale             = weibull_scale[1],
    two_side                  = two_side,
    arm2                      = arm2)

  if(arm2){
    posterior_control <- posterior_regression(
      df                        = df_c,
      family                    = family,
      alpha_max                 = alpha_max[2],
      fix_alpha                 = fix_alpha,
      prior_mean                = prior_mean,
      prior_scale               = prior_scale,
      prior_df                  = prior_df,
      prior_mean_for_intercept  = prior_mean_for_intercept,
      prior_scale_for_intercept = prior_scale_for_intercept,
      prior_df_for_intercept    = prior_df_for_intercept,
      number_mcmc               = number_mcmc,
      weibull_shape             = weibull_shape[2],
      weibull_scale             = weibull_scale[2],
      two_side                  = two_side,
      arm2                      = arm2)
  } else{
    posterior_control <- NULL
    posterior_treatment$posterior_regression$call <- formula
  }

  args1 <- list(df_t                      = df_t,
                df_c                      = NULL,
                alpha_max                 = alpha_max,
                fix_alpha                 = fix_alpha,
                prior_mean                = prior_mean,
                prior_scale               = prior_scale,
                prior_df                  = prior_df,
                prior_mean_for_intercept  = prior_mean_for_intercept,
                prior_scale_for_intercept = prior_scale_for_intercept,
                prior_df_for_intercept    = prior_df_for_intercept,
                number_mcmc               = number_mcmc,
                weibull_scale             = weibull_scale,
                weibull_shape             = weibull_shape,
                two_side                  = two_side,
                arm2                      = arm2,
                data                      = data,
                family                    = family,
                method                    = method)
  if(arm2){
    args1$df_c <- df_c
  }


  me <- list(posterior_treatment = posterior_treatment,
             posterior_control   = posterior_control,
             args1               = args1)

  class(me) <- "bdpregression"
  return(me)

})




################################################################################
# Regression posterior estimation
# 1) Estimate the discount function (if current+historical data both present)
#    - Fit two independent Bayesian regression models with "flat" priors
#    - Stochastically compare the posteriors of the historical effect
# 2) Estimate the posterior of the augmented data
#    - Augment the data by placing a prior on the current data of the mean of
#      the historical data covariate effect and the standard deviation of that
#      covariate effect, divided by alpha
################################################################################
### Combine  loss function and posterior estimation into one function
posterior_regression <- function(df, family, alpha_max, fix_alpha, prior_mean,
                                 prior_scale, prior_df, prior_mean_for_intercept,
                                 prior_scale_for_intercept, prior_df_for_intercept,
                                 number_mcmc, weibull_shape, weibull_scale,
                                 two_side, arm2){

  Y <- NULL
  historical <- NULL
  intercept <- NULL

  ### Look for "historical" column, if missing, df_ <- df and df_0 <- NULL
  hist_missing <- is.na(match("historical", colnames(df)))

  if(hist_missing){
    df_  <- df
    df_0 <- NULL
  } else{
    ### Split out historical and current data and drop historical column
    df_  <- subset(df, historical==0, select=-c(historical))  # Current data
    df_0 <- subset(df, historical==1, select=-c(historical))  # Historical data
  }

  ### If historical has a single value, then one of df_ and df_0 will have 0 rows
  if(nrow(df_) == 0){
    df_ <- NULL
  }

  if(!is.null(df_0)){
    if(nrow(df_0) == 0){
      df_0 <- NULL
    }
  }



  ##############################################################################
  # Discount function
  # - Comparison is made only if both df_ and df_0 are present
  ##############################################################################
  ### One arm discount function
  if(!arm2){
    if(is.null(df_) & is.null(df_0)){
      ### Historical and current data are both missing
      stop("Historical and/or current data was input incorrectly.")
    } else if(!is.null(df_) & !is.null(df_0)){
      ### Historical and current data are both present

      ### Extract design matrix
      X_  <- subset(df_, select = -c(Y,offset) )
      X_0 <- subset(df_0, select = -c(Y,offset) )

      ### Fit separate glm models
      posterior_flat_regression <-  bayesglm.fit(x=X_, y=df_$Y, offset=df_$offset,
        family                    = family,
        prior_mean                = prior_mean,
        prior_scale               = prior_scale,
        prior_df                  = prior_df,
        prior_mean_for_intercept  = prior_mean_for_intercept,
        prior_scale_for_intercept = prior_scale_for_intercept,
        prior_df_for_intercept    = prior_df_for_intercept)

      prior_regression <-  bayesglm.fit(x=X_0, y=df_0$Y, offset=df_0$offset,
        family                    = family,
        prior_mean                = prior_mean,
        prior_scale               = prior_scale,
        prior_df                  = prior_df,
        prior_mean_for_intercept  = prior_mean_for_intercept,
        prior_scale_for_intercept = prior_scale_for_intercept,
        prior_df_for_intercept    = prior_df_for_intercept)

      ### Set classes
      class(posterior_flat_regression) <- c("glm", "lm")
      class(prior_regression) <- c("glm", "lm")

      ### Simulate data from approximate posterior distributions of the coefficients
      sim_  <- sim(posterior_flat_regression, family=family, n.sims = number_mcmc)
      sim_0 <- sim(prior_regression, family=family, n.sims = number_mcmc)

      ### Extract the coefficient chains
      coef_  <- data.frame(sim_$coef)
      coef_0 <- data.frame(sim_0$coef)

      ### Add simulated coef chains to fit objects
      posterior_flat_regression$posterior <- coef_
      prior_regression$posterior          <- coef_0

      ### Stochastically compare the intercepts
      p_hat <- mean(coef_$intercept < coef_0$intercept)

      ### Compute mean and sd of intercepts
      mean_  <- data.frame(t(colMeans(coef_)))
      mean_0 <- data.frame(t(colMeans(coef_0)))
      sd_    <- data.frame(t(apply(coef_,2,sd)))
      sd_0   <- data.frame(t(apply(coef_0,2,sd)))

      if(fix_alpha){
        alpha_discount <- alpha_max
      } else{
        if (!two_side) {
          alpha_discount <- pweibull(p_hat, shape=weibull_shape, scale=weibull_scale)*alpha_max
        } else if (two_side){
          p_hat    <- ifelse(p_hat > 0.5, 1 - p_hat, p_hat)
          alpha_discount <- pweibull(p_hat, shape=weibull_shape, scale=weibull_scale)*alpha_max
        }
      }
    } else if(!is.null(df_) & is.null(df_0)){
      ### Historical data missing and current data present

      ### Extract design matrix
      X_  <- subset(df_, select = -c(Y,offset) )

      ### Fit separate glm models
      posterior_flat_regression <-  bayesglm.fit(x=X_, y=df_$Y, offset=df_$offset,
                                                 family                    = family,
                                                 prior_mean                = prior_mean,
                                                 prior_scale               = prior_scale,
                                                 prior_df                  = prior_df,
                                                 prior_mean_for_intercept  = prior_mean_for_intercept,
                                                 prior_scale_for_intercept = prior_scale_for_intercept,
                                                 prior_df_for_intercept    = prior_df_for_intercept)

      ### Set class
      class(posterior_flat_regression) <- c("glm", "lm")

      ### Simulate data from approximate posterior distributions of the coefficients
      sim_  <- sim(posterior_flat_regression, family=family, n.sims = number_mcmc)

      ### Extract the coefficient chains
      coef_  <- data.frame(sim_$coef)

      ### Add simulated coef chains to fit object
      posterior_flat_regression$posterior <- coef_

    } else if(is.null(df_) & !is.null(df_0)){
      ### Current data missing and historical data present

      ### Extract design matrix
      X_0 <- subset(df_0, select = -c(Y,offset) )

      ### Fit glm model
      prior_regression <-  bayesglm.fit(x=X_0, y=df_0$Y, offset=df_0$offset,
                                        family                    = family,
                                        prior_mean                = prior_mean,
                                        prior_scale               = prior_scale,
                                        prior_df                  = prior_df,
                                        prior_mean_for_intercept  = prior_mean_for_intercept,
                                        prior_scale_for_intercept = prior_scale_for_intercept,
                                        prior_df_for_intercept    = prior_df_for_intercept)

      class(prior_regression) <- c("glm", "lm")

      ### Simulate data from approximate posterior distributions of the coefficients
      sim_0 <- sim(prior_regression, family=family, n.sims = number_mcmc)

      ### Extract the coefficient chain
      coef_0 <- data.frame(sim_0$coef)

      ### Add simulated coef chains to fit object
      prior_regression$posterior          <- coef_0
    }


  } else if(arm2){   ### Two arm discount function
    if(is.null(df_) & is.null(df_0)){
      ### Historical and current data are both missing
      stop("Historical and/or current data was input incorrectly.")
    } else if(!is.null(df_) & !is.null(df_0)){
      ### Historical and current data are both present

      ### Extract design matrix
      X  <- subset(df, select=-c(historical,Y,offset))

      ### Fit separate glm models
      posterior_flat_regression <-  bayesglm.fit(x                         = X,
                                                 y                         = df$Y,
                                                 offset                    = df$offset,
                                                 family                    = family,
                                                 prior_mean                = prior_mean,
                                                 prior_scale               = prior_scale,
                                                 prior_df                  = prior_df,
                                                 prior_mean_for_intercept  = prior_mean_for_intercept,
                                                 prior_scale_for_intercept = prior_scale_for_intercept,
                                                 prior_df_for_intercept    = prior_df_for_intercept)

      class(posterior_flat_regression) <- c("glm", "lm")

      ### Simulate data from approximate posterior distributions of the coefficients
      sim_  <- sim(posterior_flat_regression, family=family, n.sims = number_mcmc)

      ### Extract the coefficient chains
      coef_  <- data.frame(sim_$coef)

      ### Add simulated coef chains to fit objects
      posterior_flat_regression$posterior <- coef_

      ### Stochastically compare the intercepts
      p_hat <- mean(coef_$treatment < 0)

      ### Compute mean and sd of covariates
      mean_  <- data.frame(t(colMeans(coef_)))
      sd_    <- data.frame(t(apply(coef_,2,sd)))

      if(fix_alpha){
        alpha_discount <- alpha_max
      } else{
        if (!two_side) {
          alpha_discount <- pweibull(p_hat, shape=weibull_shape, scale=weibull_scale)*alpha_max
        } else if (two_side){
          p_hat    <- ifelse(p_hat > 0.5, 1 - p_hat, p_hat)
          alpha_discount <- pweibull(p_hat, shape=weibull_shape, scale=weibull_scale)*alpha_max
        }
      }
    }


  }


  ##############################################################################
  # Posterior augmentation
  # - If current or historical data are missing, this will not augment
  ##############################################################################
  ### One arm data augmentation
  ### - Augment via a prior on the intercept
  if(!arm2){
    if(!is.null(df_) & !is.null(df_0)){

      # Parse out covariate estimates
      prior_mean  <- as.numeric(subset(mean_0, select = -intercept))
      prior_scale <- as.numeric(subset(sd_0, select = -intercept))

      posterior_regression <-  bayesglm.fit(x=X_, y=df_$Y, offset=df_$offset,
        family                    = family,
        prior_mean                = prior_mean,
        prior_scale               = prior_scale/sqrt(alpha_discount),
        prior_df                  = prior_df,
        prior_mean_for_intercept  = mean_0$intercept,
        prior_scale_for_intercept = sd_0$intercept/sqrt(alpha_discount),
        prior_df_for_intercept    = prior_df_for_intercept)


      class(posterior_regression) <- c("glm", "lm")

      ### Simulate date from approximate posterior distribution of the coefficients
      sim_a  <- sim(posterior_regression, family=family, n.sims = number_mcmc)

      ### Extract the coefficient chain
      coef_a  <- data.frame(sim_a$coef)

      ### Add simulated coef chains to fit object
      posterior_regression$posterior <- coef_a

      ### Collect output to save, assign new class and return
      res <- list(p_hat                     = p_hat,
                  alpha_discount            = alpha_discount,
                  posterior_regression      = posterior_regression,
                  posterior_flat_regression = posterior_flat_regression,
                  prior_regression          = prior_regression)

      return(res)
    } else if(!is.null(df_) & is.null(df_0)){

      res <- list(p_hat                     = NULL,
                  alpha_discount            = NULL,
                  posterior_regression      = posterior_flat_regression,
                  posterior_flat_regression = posterior_flat_regression,
                  prior_regression          = NULL)
      return(res)

    } else if(is.null(df_) & !is.null(df_0)){

      res <- list(p_hat                     = NULL,
                  alpha_discount            = NULL,
                  posterior_regression      = prior_regression,
                  posterior_flat_regression = NULL,
                  prior_regression          = prior_regression)
      return(res)
    }


  } else if(arm2){


  }

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




bayesglm.fit <- function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
                          mustart = NULL, offset = rep(0, nobs), family = gaussian(),
                          control = list(), intercept = TRUE,
                          prior_mean = 0, prior_scale = NULL, prior_df = 1,
                          prior_mean_for_intercept = 0,
                          prior_scale_for_intercept = NULL,
                          prior_df_for_intercept = 1,
                          min_prior_scale = 1e-12, scaled = TRUE,
                          print_unnormalized_log_posterior = FALSE,
                          Warning = TRUE){

  control <- do.call("glm.control", control)

  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y)){
    rownames(y)
  }else{
    names(y)
  }
  conv <- FALSE
  nobs <- NROW(y)
  nvars <- NCOL(x)
  #===============================
  #  initialize priors
  #===============================
  if(is.null(prior_scale)){
    prior_scale <- 2.5
    if(family$link == "probit"){
      prior_scale <- prior_scale*1.6
    }
  }

  if(is.null(prior_scale_for_intercept)){
    prior_scale_for_intercept <- 10
    if(family$link == "probit"){
      prior_scale_for_intercept <- prior_scale_for_intercept*1.6
    }
  }

  if(intercept){
    nvars <- nvars - 1
  }

  if(length(prior_mean)==1L){
    prior_mean <- rep(prior_mean, nvars)
  }else if(length(prior_mean)!=nvars){
    stop("invalid length for prior_mean")
  }

  if(length(prior_scale)==1L){
    prior_scale <- rep(prior_scale, nvars)
  }else if(length(prior_scale)!=nvars){
    stop("invalid length for prior_scale")
  }

  if(length(prior_df)==1L){
    prior_df <- rep(prior_df, nvars)
  }else if(length(prior_df)!=nvars){
    stop("invalid length for prior_df")
  }

  if(intercept){
    prior_mean <- c(prior_mean_for_intercept, prior_mean)
    prior_scale <- c(prior_scale_for_intercept, prior_scale)
    prior_df <- c(prior_df_for_intercept, prior_df)
  }

  if(scaled){
    if(family$family=="gaussian"){
      prior_scale <- prior_scale*2*sd(y)
    }
    prior_scale_0 <- prior_scale
    if(nvars==0) nvars = 1
    for(j in 1:nvars){
      x.obs <- x[,j]
      x.obs <- x.obs[!is.na(x.obs)]
      num.categories <- length(unique(x.obs))
      x.scale <- 1
      if(num.categories==2L){
        x.scale <- max(x.obs) - min(x.obs)
      }else if(num.categories>2){
        x.scale <- 2*sd(x.obs)
      }
      prior_scale[j] <- prior_scale[j]/x.scale
      if(prior_scale[j] < min_prior_scale){
        prior_scale[j] <- min_prior_scale
        warning("prior scale for varible ", j,
                " set to min_prior_scale = ", min_prior_scale, "\n")
      }
    }
  }

  #===================
  nvars <- NCOL(x)
  EMPTY <- nvars == 0
  if (is.null(weights))
    weights <- rep.int(1, nobs)
  if (is.null(offset))
    offset <- rep.int(0, nobs)
  variance <- family$variance
  linkinv <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv))
    stop("'family' argument seems not to be a valid family object",
         call. = FALSE)
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  unless.null <- function(x, if.null){
    if (is.null(x))
      if.null
    else x
  }
  valideta <- unless.null(family$valideta, function(eta) TRUE)
  validmu <- unless.null(family$validmu, function(mu) TRUE)
  if (is.null(mustart)) {
    eval(family$initialize)
  }else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  if (EMPTY) {
    eta <- rep.int(0, nobs) + offset
    if (!valideta(eta))
      stop("invalid linear predictor values in empty model",
           call. = FALSE)
    mu <- linkinv(eta)
    if (!validmu(mu))
      stop("invalid fitted means in empty model", call. = FALSE)
    dev <- sum(dev.resids(y, mu, weights))
    w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
    residuals <- (y - mu)/mu.eta(eta)
    good <- rep_len(TRUE, length(residuals))
    boundary <- conv <- TRUE
    coef <- numeric()
    iter <- 0L
  } else {
    coefold <- NULL
    eta <- if (!is.null(etastart)){
      etastart
    }else if (!is.null(start)){
      if (length(start) != nvars){
        if(start==0&length(start)==1){
          start <- rep(0, nvars)
          offset + as.vector(ifelse((NCOL(x) == 1L), x*start, x %*% start))
        }else{
          stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                        nvars, paste(deparse(xnames), collapse = ", ")),
               domain = NA)
        }
      } else {
        coefold <- start
        offset + as.vector(if (NCOL(x) == 1L)
          x * start
          else x %*% start)
      }
    }else{
      family$linkfun(mustart)
    }
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta)))
      stop("cannot find valid starting values: please specify some",
           call. = FALSE)
    devold <- sum(dev.resids(y, mu, weights))
    boundary <- conv <- FALSE

    #======================================
    #   initialize prior_sd
    #======================================
    prior_sd <- prior_scale
    #=====================================
    dispersion <- ifelse((family$family %in% c("poisson", "binomial")),  1, var(y)/10000)
    dispersionold <- dispersion
    for (iter in 1L:control$maxit) {
      good <- weights > 0
      varmu <- variance(mu)[good]
      if (anyNA(varmu))
        stop("NAs in V(mu)")
      if (any(varmu == 0))
        stop("0s in V(mu)")
      mu.eta.val <- mu.eta(eta)
      if (any(is.na(mu.eta.val[good])))
        stop("NAs in d(mu)/d(eta)")
      good <- (weights > 0) & (mu.eta.val != 0)
      if (all(!good)) {
        conv <- FALSE
        warning(gettextf("no observations informative at iteration %d",
                         iter), domain = NA)
        break
      }
      z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
      w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
      ngoodobs <- as.integer(nobs - sum(!good))
      #======================
      #  data augmentation
      #=========================
      x.star <- rbind(x, diag(NCOL(x)))
      if(intercept&scaled){
        x.star[nobs+1,] <- colMeans(x)
      }
      z.star <- c (z, prior_mean)
      w.star <- c (w, sqrt(dispersion)/prior_sd)
      #=================================================
      good.star <- c (good, rep(TRUE,NCOL(x)))
      ngoodobs.star <- ngoodobs + NCOL(x)
      fit <- lm.fit(x = x.star[good.star,,drop=FALSE]*w.star, y = z.star*w.star)
      if (any(!is.finite(fit$coefficients))) {
        conv <- FALSE
        warning(gettextf("non-finite coefficients at iteration %d",
                         iter), domain = NA)
        break
      }
      start[fit$qr$pivot] <- coefs.hat <- fit$coefficients
      fit$qr$qr <- as.matrix (fit$qr$qr)
      V.coefs <- chol2inv(fit$qr$qr[1:NCOL(x.star), 1:NCOL(x.star), drop = FALSE])
      if (family$family == "gaussian" & scaled){
        prior_scale <- prior_scale_0
      }
      prior_sd <- ifelse(prior_df == Inf, prior_scale,
                         sqrt(((coefs.hat - prior_mean)^2 + diag(V.coefs)*dispersion +
                                 prior_df * prior_scale^2)/(1 + prior_df)))
      start[fit$qr$pivot] <- fit$coefficients
      eta <- drop(x %*% start)
      mu <- linkinv(eta <- eta + offset)
      dev <- sum(dev.resids(y, mu, weights))
      if (!(family$family %in% c("poisson", "binomial"))) {
        mse.resid <- mean((w * (z - x %*% coefs.hat))^2)
        mse.uncertainty <- mean(rowSums(( x %*% V.coefs ) * x)) * dispersion # faster
        dispersion <- mse.resid + mse.uncertainty
      }
      if (control$trace)
        cat("Deviance = ", dev, " Iterations - ", iter,
            "\n", sep = "")
      boundary <- FALSE
      if (!is.finite(dev)) {
        if (is.null(coefold))
          stop("no valid set of coefficients has been found: please supply starting values",
               call. = FALSE)
        warning("step size truncated due to divergence",
                call. = FALSE)
        ii <- 1
        while (!is.finite(dev)) {
          if (ii > control$maxit)
            stop("inner loop 1; cannot correct step size",
                 call. = FALSE)
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
          dev <- sum(dev.resids(y, mu, weights))
        }
        boundary <- TRUE
        if (control$trace)
          cat("Step halved: new deviance = ", dev, "\n",
              sep = "")
      }
      if (!(valideta(eta) && validmu(mu))) {
        if (is.null(coefold))
          stop("no valid set of coefficients has been found: please supply starting values",
               call. = FALSE)
        warning("step size truncated: out of bounds",
                call. = FALSE)
        ii <- 1
        while (!(valideta(eta) && validmu(mu))) {
          if (ii > control$maxit)
            stop("inner loop 2; cannot correct step size",
                 call. = FALSE)
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
        }
        boundary <- TRUE
        dev <- sum(dev.resids(y, mu, weights))
        if (control$trace)
          cat("Step halved: new deviance = ", dev, "\n",
              sep = "")
      }
      #===============================
      # print unnormalized log posterior
      #================================
      if (family$family == "binomial" && print_unnormalized_log_posterior) {
        logprior <- sum(dt(coefs.hat, prior_df, prior_mean, log = TRUE))
        xb <- 1/(1 + exp(-( x %*% coefs.hat )))
        loglikelihood <- sum( log( c( xb[ y == 1 ], 1 - xb[ y == 0 ] ) ) )
        cat( "log prior: ", logprior, ", log likelihood: ", loglikelihood, ",
             unnormalized log posterior: ", loglikelihood +logprior, "\n" ,sep="")
      }
      #================================

      if (iter > 1 & abs(dev - devold)/(0.1 + abs(dev)) <
          control$epsilon & abs(dispersion - dispersionold)/(0.1 +
                                                             abs(dispersion)) < control$epsilon) {
        conv <- TRUE
        coef <- start
        break
      }else {
        devold <- dev
        dispersionold <- dispersion
        coef <- coefold <- start
      }
      }
    if (!conv){
      warning("algorithm did not converge", call. = FALSE)
    }
    if (boundary){
      warning("algorithm stopped at boundary value",
              call. = FALSE)
    }
    eps <- 10 * .Machine$double.eps
    if (family$family == "binomial") {
      if (any(mu > 1 - eps) || any(mu < eps)){
        warning("fitted probabilities numerically 0 or 1 occurred",
                call. = FALSE)
      }
    }
    if (family$family == "poisson") {
      if (any(mu < eps)){
        warning("fitted rates numerically 0 occurred",
                call. = FALSE)
      }
    }
    if (fit$rank < nvars){
      coef[fit$qr$pivot][seq.int(fit$rank + 1, nvars)] <- NA
    }
    xxnames <- xnames[fit$qr$pivot]
    residuals <- rep.int(NA, nobs)
    residuals[good] <- z - (eta - offset)[good]
    fit$qr$qr <- as.matrix(fit$qr$qr)
    nr <- min(sum(good), nvars)
    if (nr < nvars) {
      Rmat <- diag(nvars)
      Rmat[1L:nr, 1L:nvars] <- fit$qr$qr[1L:nr, 1L:nvars]
    } else Rmat <- fit$qr$qr[1L:nvars, 1L:nvars]
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    names(coef) <- xnames
    colnames(fit$qr$qr) <- xxnames
    dimnames(Rmat) <- list(xxnames, xxnames)
  }
  names(residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames
  wt <- rep.int(0, nobs)
  wt[good] <- w^2
  names(wt) <- ynames
  names(weights) <- ynames
  names(y) <- ynames
  wtdmu <- if (intercept){
    sum(weights * y)/sum(weights)
  } else{
    linkinv(offset)
  }
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)
  rank <- if (EMPTY) {
    0
  } else{
    fit$rank
  }
  resdf <- n.ok - rank
  aic.model <- aic(y, n.ok, mu, weights, dev) + 2 * rank
  list(coefficients = coef,
       residuals = residuals,
       fitted.values = mu,
       effects = if (!EMPTY) fit$effects,
       R = if (!EMPTY) Rmat,
       rank = rank,
       qr = if (!EMPTY) structure(getQr(fit)[c("qr", "rank", "qraux", "pivot", "tol")], class = "qr"),
       family = family,
       linear.predictors = eta,
       deviance = dev,
       aic = aic.model,
       null.deviance = nulldev,
       iter = iter,
       weights = wt,
       prior.weights = weights,
       df.residual = resdf,
       df.null = nulldf,
       y = y,
       converged = conv,
       boundary = boundary,
       prior_mean = prior_mean,
       prior_scale = prior_scale,
       prior_df = prior_df,
       prior_sd = prior_sd,
       dispersion = dispersion)
}

getQr <- function(x, ...){
  if (is.null(r <- x$qr))
    stop("lm object does not have a proper 'qr' component.\n Rank zero or should not have used lm(.., qr=FALSE).")
  r
}


mvrnorm <- function (n = 1, mu, Sigma, tol = 1e-06){
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p)))
      stop("incompatible arguments")
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L])))
      stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*%
      t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
      nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1)
      drop(X)
  else t(X)
}



sim <- function(object, family=gaussian, n.sims=100){

  if(family$family == "gaussian"){
    coef      <- object$coef
    sigma.hat <- object$dispersion
    beta.hat  <- coef
    V.beta    <- summary(object)$cov.unscaled
    n         <- length(object$y)
    k         <- object$rank

    sigma <- rep (NA, n.sims)
    beta  <- array (NA, c(n.sims,k))
    dimnames(beta) <- list(NULL, names(beta.hat))

    for (s in 1:n.sims){
      sigma[s] <- sigma.hat*sqrt((n-k)/rchisq(1,n-k))
      beta[s,] <- mvrnorm(1, beta.hat, V.beta*sigma[s]^2)
    }

    ans <- list(coef = beta, sigma = sigma)
    return(ans)
  } else{
    summ      <- summary(object, correlation=TRUE, dispersion=object$dispersion)
    coef      <- summ$coef
    beta.hat  <- coef[,1,drop=FALSE]
    sd.beta   <- coef[,2,drop=FALSE]
    corr.beta <- summ$corr
    n         <- summ$df[1] + summ$df[2]
    k         <- summ$df[1]
    V.beta    <- corr.beta*array(sd.beta,c(k,k))*t(array(sd.beta,c(k,k)))
    beta      <- array (NA, c(n.sims,k))
    dimnames(beta) <- list (NULL, dimnames(beta.hat)[[1]])

    for (s in 1:n.sims){
      beta[s,] <- mvrnorm(1, beta.hat, V.beta)
    }

    beta2           <- array (0, c(n.sims,length(coefficients(object))))
    dimnames(beta2) <- list (NULL, names(coefficients(object)))
    beta2[,dimnames(beta2)[[2]]%in%dimnames(beta)[[2]]] <- beta

    sigma <- rep(sqrt(summ$dispersion), n.sims)
    ans   <- list(coef = beta2, sigma = sigma)
    return(ans)
  }
}
