## ---- SETTINGS-knitr, include=FALSE--------------------------------------
library(bayesDP)
stopifnot(require(knitr))
opts_chunk$set(
  #comment = NA,
  #message = FALSE,
  #warning = FALSE,
  #eval = if (isTRUE(exists("params"))) params$EVAL else FALSE,
  dev = "png",
  dpi = 150,
  fig.asp = 0.8,
  fig.width = 5,
  out.width = "60%",
  fig.align = "center"
  )

# Run two models to document the discount function plots
fit01 <- bdpnormal(mu_t=10, sigma_t = 1, N_t=500, 
                   mu0_t=10, sigma0_t = 1, N0_t=500)
fit02 <- bdpnormal(mu_t=10, sigma_t = 1, N_t=500, 
                   mu0_t=10, sigma0_t = 1, N0_t=500)

fit_scaledweibull <- bdpbinomial(y_t=10, N_t=500, y0_t=25, N0_t=250, 
                                 discount_function="scaledweibull")
fit_identity <- bdpbinomial(y_t=10, N_t=500, y0_t=25, N0_t=250,
                            discount_function="identity")

## ---- echo=FALSE---------------------------------------------------------
df1 <- plot(fit02, type="discount", print=FALSE)
df1 + ggtitle("Discount function plot", "Weibull distribution with shape=3 and scale=0.135")

## ---- echo=FALSE---------------------------------------------------------
df2 <- plot(fit_identity, type="discount", print=FALSE)
df2 + ggtitle("Discount function plot", "Identity")

## ------------------------------------------------------------------------
p1 <- plot(fit01, type="discount", print=FALSE)
p1 + ggtitle("Discount Function Plot :-)")

## ------------------------------------------------------------------------
set.seed(42)
fit1 <- bdpnormal(mu_t      = 30,
                  sigma_t   = 10,
                  N_t       = 250, 
                  mu0_t     = 50,
                  sigma0_t  = 5,
                  N0_t      = 250,
                  alpha_max = 1,
                  fix_alpha = TRUE)
summary(fit1)

## ------------------------------------------------------------------------
set.seed(42)
fit1a <- bdpnormal(mu_t     = 30,
                   sigma_t   = 10,
                   N_t       = 250, 
                   mu0_t     = 50,
                   sigma0_t  = 5,
                   N0_t      = 250,
                   fix_alpha = FALSE)
summary(fit1a)

## ------------------------------------------------------------------------
mean_augmented <- round(median(fit1a$posterior_treatment$posterior_mu),4)
CI95_augmented <- round(quantile(fit1a$posterior_treatment$posterior_mu, prob=c(0.025, 0.975)),4)

## ------------------------------------------------------------------------
plot(fit1a)

## ------------------------------------------------------------------------
set.seed(42)
fit2 <- bdpnormal(mu_t      = 30,
                  sigma_t   = 10,
                  N_t       = 250, 
                  mu0_t     = 50,
                  sigma0_t  = 5,
                  N0_t      = 250,
                  mu_c      = 25,
                  sigma_c   = 10,
                  N_c       = 250, 
                  mu0_c     = 25,
                  sigma0_c  = 5,
                  N0_c      = 250,
                  fix_alpha = FALSE)
summary(fit2)

## ------------------------------------------------------------------------
plot(fit2)

