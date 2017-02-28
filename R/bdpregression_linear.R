#' bdpregression_linear
#' bdpregression_linear
#' @title bdpregression_linear: Post Estimate Gaussian
#' @param data numeric number
#' @param formula formula
#' @param family family
#' @param treatment name of the treatment effect in the data dataframe.
#' @param mu0 prior mean of the treatment effect
#' @param sigma02 prior variance of the treatment effect.
#' @param prior.dist prior distribution of the treatment and covariate effects.
#' @param prior.mean prior mean of the covariate effects. Default is zero.
#' @param prior.scale prior variance of the covariate effects. Default is 1e6.
#' @param prior.df prior degrees of freedom of the treatment and covariate effects.
#' @param weibull_scale weibull_scale
#' @param weibull_shape weibull_shape
#' @param alpha_max alpha_max
#' @param number_mcmc number_mcmc
#' @param two_side numeric
#' @examples
#' set.seed(42)
#' data <- data.frame(y         = rnorm(100, 4, 0.1),
#'                    x         = c(rnorm(50,1,0.1), rnorm(50,3,0.1)),
#'                    treatment = c(rep(0,50),rep(1,50)))
#'
#' fit <- bdpregression_linear(data,
#' formula       = y ~ treatment + x,
#' family        = "gaussian",
#' treatment     = "treatment",
#' prior.dist    = NULL,
#' prior.mean    = 0,
#' prior.scale   = 1000,
#' prior.df      = Inf,
#' mu0           = 1,
#' sigma02       = 0.1,
#' weibull_scale = 1,
#' weibull_shape = 1,
#' alpha_max     = 1,
#' number_mcmc   = 10000,
#' two_side      = 0)
#' ### Main parameter of interest:
#' fit$effect_est
#' @rdname bdpregression_linear
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @aliases bdpregression_linear,ANY-method
#' @export bdpregression_linear
bdpregression_linear <- setClass("bdpregression_linear", slots = c(est = "list",
                                                                   f1 = "list",
                                                                   args1 = "list"))

############################################################################
### Need to fix plots: should not involve effective sample sizes, instead,
### weight of alpha from 0 --> 1

setGeneric("bdpregression_linear",
           function(data,
                    formula       = y ~ treatment + x,
                    family        = "gaussian",
                    treatment     = "treatment",
                    prior.dist    = NULL,
                    prior.mean    = 0,
                    prior.scale   = 1000,
                    prior.df      = Inf,
                    mu0           = 1,
                    sigma02       = 0.1,
                    weibull_scale = 1,
                    weibull_shape = 1,
                    alpha_max     = 1,
                    number_mcmc   = 10000,
                    two_side      = 0){
             standardGeneric("bdpregression_linear")
           })

setMethod("bdpregression_linear",
          signature(),
          function(data,
                   formula       = y ~ treatment + x,
                   family        = "gaussian",
                   treatment     = "treatment",
                   prior.dist    = NULL,
                   prior.mean    = 0,
                   prior.scale   = 1000,
                   prior.df      = Inf,
                   mu0           = 1,
                   sigma02       = 0.1,
                   weibull_scale = 1,
                   weibull_shape = 1,
                   alpha_max     = 1,
                   number_mcmc   = 10000,
                   two_side      = 0){

  N_new        <- nrow(data)


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


  ################################################################################
  # Estimate alpha, discount function
  # - Currently, only place prior on (adjusted) treatment effect
  ################################################################################
  Loss_function <- function(Y, X, offset, family, treatment, mu0, sigma02, alpha_max,
                            prior.df, prior.scale, prior.mean,
                            number_mcmc, weibull_shape, weibull_scale, two_side){


    ### Get treatment effect estimate using flat prior
    q_flat <- bayesglm.fit(x           = X,
                           y           = Y,
                           offset      = offset,
                           prior.df    = prior.df,
                           prior.scale = sqrt(prior.scale),
                           prior.mean  = prior.mean,
                           family      = family)
    class(q_flat) <- c("bayesglm", "glm", "lm")


    ### Grab mean/sd of treatment effect
    coef_table <- summary(q_flat)$coefficients
    coef_trt   <- which(row.names(coef_table) == treatment)
    mu_new     <- coef_table[coef_trt,1]
    sd_new     <- coef_table[coef_trt,2]


    ### Compare new data to "historical data"
    mu_post_flat  <- rnorm(number_mcmc, mu_new, sd_new)
    mu_post_flat0 <- rnorm(number_mcmc, mu0, sigma02)

    ### Test of model vs real
    p_test  <- mean(mu_post_flat < mu_post_flat0)

    ### Number of effective sample size given shape and scale loss function
    if (two_side == 0) {
      alpha_loss <- pweibull(p_test, shape = weibull_shape, scale = weibull_scale) * alpha_max
    } else if (two_side == 1){
      p_test1    <- ifelse(p_test > 0.5, 1 - p_test, p_test)
      alpha_loss <- pweibull(p_test1, shape = weibull_shape, scale = weibull_scale) * alpha_max
    }

    return(list(alpha_loss   = alpha_loss,
                pvalue       = p_test,
                mu_post_flat = mu_post_flat,
                mu0          = mu_post_flat0))
  }


  ### Scale the standard deviation of the historical group by alpha0
  #std0        <- prior_sigma/sqrt(alpha0)


  ################################################################################
  # Estimate posterior for (adjusted) treatment effect given alpha_loss value
  ################################################################################
  mu_post_aug <- function(Y, X, offset, family, treatment, mu0, sigma02, alpha_loss,
                          prior.df, prior.scale, prior.mean,
                          number_mcmc, weibull_shape, weibull_scale) {

    ### Scale the historicall variance with respect to the alpha value
    std0        <- sqrt(sigma02)/sqrt(alpha_loss)

    ### Match historical treatment info to the correct slot
    nP                 <- ncol(X)-1
    namesX             <- colnames(X)[-1]
    idTrt              <- which(namesX == treatment)
    prior_mean         <- rep(prior.mean, nP)
    prior_mean[idTrt]  <- mu0
    prior_scale        <- rep(prior.scale, nP)
    prior_scale[idTrt] <- std0^2


    q_final  <- bayesglm.fit(x           = X,
                             y           = Y,
                             offset      = offset,
                             prior.df    = prior.df,
                             prior.scale = sqrt(prior_scale),
                             prior.mean  = prior_mean,
                             family      = family)
    class(q_final) <- c("bayesglm", "glm", "lm")


    ### Grab mean/sd of treatment effect
    coef_table <- summary(q_final)$coefficients
    coef_trt   <- which(row.names(coef_table) == treatment)
    mu_new     <- coef_table[coef_trt,1]
    sd_new     <- coef_table[coef_trt,2]


    ### Calculate posterior of the treatment effect
    mu_post <- rnorm(number_mcmc, mu_new, sd_new)


    ### Extract posterior estimates of intercept, treatment effect, and covariate(s)
    mu_covariate_post <- summary(q_final)$coefficients[,1]
    sd_covariate_post <- vcov(q_final)


    ### Draw posterior samples of the regression estimates
    ### - This is a Laplace approximation - we will derive one other option
    covariate_post <- MASS::mvrnorm(number_mcmc, mu_covariate_post, sd_covariate_post)


    return(list(mu_post           = mu_post,
                mu_covariate_post = mu_covariate_post,
                sd_covariate_post = sd_covariate_post,
                covariate_post    = covariate_post))
  }


  ################################################################################
  # Combine loss function and posterior estimation into one function             #
  ################################################################################
  mu_posterior <- function(Y, X, offset, family, treatment, mu0, sigma02, alpha_max,
                           prior.df, prior.scale, prior.mean,
                           number_mcmc, weibull_shape, weibull_scale, two_side) {

    alpha_loss <- Loss_function(Y             = Y,
                                X             = X,
                                offset        = offset,
                                family        = family,
                                treatment     = treatment,
                                mu0           = mu0,
                                sigma02       = sigma02,
                                alpha_max     = alpha_max,
                                prior.df      = prior.df,
                                prior.scale   = prior.scale,
                                prior.mean    = prior.mean,
                                number_mcmc   = number_mcmc,
                                weibull_shape = weibull_shape,
                                weibull_scale = weibull_scale,
                                two_side      = two_side)

    mu_posterior <- mu_post_aug(Y             = Y,
                                X             = X,
                                offset        = offset,
                                family        = family,
                                treatment     = treatment,
                                mu0           = mu0,
                                sigma02       = sigma02,
                                alpha_loss    = alpha_loss$alpha_loss,
                                prior.df      = prior.df,
                                prior.scale   = prior.scale,
                                prior.mean    = prior.mean,
                                number_mcmc   = number_mcmc,
                                weibull_shape = weibull_shape,
                                weibull_scale = weibull_scale)


    return(list(alpha_loss        = alpha_loss$alpha_loss,
                pvalue            = alpha_loss$pvalue,
                mu_posterior      = mu_posterior$mu_post,
                mu_posterior_flat = alpha_loss$mu_post_flat,
                mu_prior          = alpha_loss$mu0,
                weibull_scale     = weibull_scale,
                weibull_shape     = weibull_shape,
                mu0               = mu0,
                sigma02           = sigma02))
  }


  final <- function(posterior) {
    den_post  <- density(posterior$mu_posterior, adjust = 0.5)
    den_flat  <- density(posterior$mu_posterior_flat, adjust = 0.5)
    den_prior <- density(posterior$mu_prior, adjust = 0.5)

    post <- posterior$mu_posterior

    return(list(den_post  = den_post,
                den_flat  = den_flat,
                den_prior = den_prior,
                post      = post))
  }


  est <- mu_posterior(
    Y             = Y,
    X             = X,
    offset        = offset,
    family        = family,
    treatment     = treatment,
    mu0           = mu0,
    sigma02       = sigma02,
    alpha_max     = alpha_max,
    prior.df      = prior.df,
    prior.scale   = prior.scale,
    prior.mean    = prior.mean,
    number_mcmc   = number_mcmc,
    weibull_shape = weibull_shape,
    weibull_scale = weibull_scale,
    two_side      = two_side)            #Two or one sided hypothesis test?


  f1 <- final(posterior = est)

  args1 <- list(data          = data,
                formula       = formula,
                family        = family,
                treatment     = treatment,
                prior.dist    = prior.dist,
                prior.mean    = prior.mean,
                prior.scale   = prior.scale,
                prior.df      = prior.df,
                mu0           = mu0,
                sigma02       = sigma02,
                weibull_scale = weibull_scale,
                weibull_shape = weibull_shape,
                alpha_max     = alpha_max,
                number_mcmc   = number_mcmc,
                two_side      = two_side)

  me <- list(est   = est,
             f1    = f1,
             args1 = args1)

  class(me) <- "bdpregression_linear"

  return(me)

})
