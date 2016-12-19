## declare the display generic
setGeneric("post_est_gaussian", function(data,
                                         formula,
                                         prior_mu,
                                         prior_sigma,
                                         weibull_scale,
                                         weibull_shape,
                                         alpha_max,
                                         N0_max = NULL,
                                         N_mcmc)
  standardGeneric("post_est_gaussian")
)

#setMethod("post_est_gaussian",
#          signature(data = "ANY"),
#          function(data){
#            message("Wrong object")
#          })

#setMethod("post_est_gaussian",
#          signature(data = "missing"),
#          function(data){
#            message("Missing object")
#          })

#' Post Estimate Gaussian
#'
#' Post Estimate Gaussian
#'
#' @title post_est_gaussian: Post Estimate Gaussian
#' @param data numeric number
#' @param formula formula
#' @param prior_mu prior_mu
#' @param prior_sigma prior_sigma
#' @param weibull_scale weibull_scale
#' @param weibull_shape weibull_shape
#' @param alpha_max alpha_max
#' @param N0_max N0_max
#' @param N_mcmc N_mcmc
#'
#' @examples
#' set.seed(42)
#' data <- data.frame(y         = rnorm(100, 4, 0.1),
#'                    x         = c(rnorm(50,1,0.1), rnorm(50,3,0.1)),
#'                    treatment = c(rep(0,50),rep(1,50)))
#'
#' fit <- post_est_gaussian(data          = data,
#'                          formula       = y ~ treatment + x,
#'                          prior_mu      = 1,
#'                          prior_sigma   = 0.1,
#'                          weibull_scale = 1,
#'                          weibull_shape = 1,
#'                          alpha_max     = 1,
#'                          N_mcmc        = 10000)
#'
#' ### Main parameter of interest:
#' fit$effect_est
#'
#' @rdname post_est_gaussian
#' @export post_est_gaussian

setMethod("post_est_gaussian",
          signature(data = "data.frame",
                    formula = "formula",
                    prior_mu = "numeric",
                    prior_sigma = "numeric",
                    weibull_scale = "numeric",
                    weibull_shape = "numeric",
                    alpha_max = "numeric",
                    N0_max = "missing",
                    N_mcmc = "numeric"),

          function(data,
                   formula,
                   prior_mu,
                   prior_sigma,
                   weibull_scale,
                   weibull_shape,
                   alpha_max,
                   N0_max = NULL,
                   N_mcmc){

N_new        <- nrow(data)

################################################################################
### Estimate alpha (loss function)
################################################################################
q_flat       <- bayesglm(formula,
                         prior.df    = Inf,
                         prior.scale = 1000,
                         prior.mean  = 0,
                         family      = gaussian,
                         data        = data)

### Grab mean/sd of treatment effect
mu_new      <- summary(q_flat)$coefficients[2,1]
sd_new      <- summary(q_flat)$coefficients[2,2]

### Compare new data to "historical data"
mc_new_data  <- rnorm(N_mcmc,mu_new,sd_new)
mc_hist_data <- rnorm(N_mcmc, prior_mu, prior_sigma)

### Compare distributions of current and historical data
p           <- mean(mc_hist_data<mc_new_data)
p_test1     <- ifelse(p > 0.5, 1 - p, p)

### Compute alpha0, loss function, from Weibull loss function
alpha0      <- pweibull(p_test1, shape = weibull_shape, scale = weibull_scale)*alpha_max
alpha0      <- ifelse(alpha0<1e-12,1e-12,alpha0)

### Scale the standard deviation of the historical group by alpha0
std0        <- prior_sigma/sqrt(alpha0)

### Note: if N0_max is not null, then use
if(!is.null(N0_max)){
  alpha0     <- pweibull(p_test1, shape = weibull_shape, scale = weibull_scale)
  alpha0     <- ifelse(alpha0<1e-12,1e-12,alpha0)
  N_eff      <- N0_max*alpha0
  std0       <- prior_sigma/sqrt(N0_max)
}

################################################################################
### Estimate treatment effect of current data with historical prior,
### where historical prior has been down-weighted using alpha0
################################################################################
q_final     <- bayesglm(formula,
                        prior.df    = Inf,
                        prior.scale = c(std0, 1000),
                        prior.mean  = c(prior_mu, 0),
                        family      = gaussian,
                        data        = data)

### Extract posterior estimates of intercept, treatment effect, and covariate x
mu_post       <- summary(q_final)$coefficients[,1]
sd_post       <- vcov(q_final)

### Draw posterior samples of the regression estimates
### - This is a Laplace approximation - we will derive one other option
mc_final_data <- mvrnorm(N_mcmc, mu_post, sd_post)

out <- list(alpha0           = alpha0,
            bayesian_p_value = p,
            N_new            = N_new,
            mu_post          = mu_post,
            sd_post          = sd_post,
            mc_final_data    = mc_final_data,
            mc_hist_data     = mc_hist_data,
            mc_new_data      = mc_new_data,
            effect_est       = coef(q_final)[2])

return(out)})
