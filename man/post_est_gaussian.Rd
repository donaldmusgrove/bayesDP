set.seed(42)
data <- data.frame(y         = rnorm(100, 4, 0.1),
                   x         = c(rnorm(50,1,0.1), rnorm(50,3,0.1)),
                   treatment = c(rep(0,50),rep(1,50)))

fit <- post_est_gaussian(data          = data,
                         prior_mu      = 1,
                         prior_sigma   = 0.1,
                         weibull_scale = 1,
                         weibull_shape = 1,
                         alpha_max     = 1,
                         N_mcmc        = 10000)

### Main parameter of interest:
fit$effect_est
