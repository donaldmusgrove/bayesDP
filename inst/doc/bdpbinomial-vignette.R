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
fit01 <- bdpbinomial(y_t=10, N_t=500, y0_t=25, N0_t=250, two_side=FALSE)
fit02 <- bdpbinomial(y_t=10, N_t=500, y0_t=25, N0_t=250)


## ---- echo=FALSE---------------------------------------------------------
plot(fit02, type="discount")

## ---- echo=FALSE---------------------------------------------------------
plot(fit01, type="discount")

## ------------------------------------------------------------------------
set.seed(42)
fit1 <- bdpbinomial(y_t       = 10,
                    N_t       = 200,
                    y0_t      = 25,
                    N0_t      = 250,
                    alpha_max = 1,
                    fix_alpha = TRUE)
summary(fit1)

## ------------------------------------------------------------------------
set.seed(42)
fit1a <- bdpbinomial(y_t       = 10,
                     N_t       = 200,
                     y0_t      = 25,
                     N0_t      = 250,
                     alpha_max = 1,
                     fix_alpha = FALSE)
summary(fit1a)

## ------------------------------------------------------------------------
mean_augmented <- round(median(fit1a$posterior_treatment$posterior),4)
CI95_augmented <- round(quantile(fit1a$posterior_treatment$posterior, prob=c(0.025, 0.975)),4)

## ------------------------------------------------------------------------
plot(fit1a)

## ------------------------------------------------------------------------
set.seed(42)
fit2 <- bdpbinomial(y_t  = 10,
                    N_t  = 200,
                    y0_t = 25,
                    N0_t = 250,
                    y_c  = 15,
                    N_c  = 200,
                    y0_c = 20,
                    N0_c = 250)
summary(fit2)

## ------------------------------------------------------------------------
plot(fit2)

