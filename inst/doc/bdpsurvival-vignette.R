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
time   <- c(rexp(50, rate=1/20), rexp(50, rate=1/10))
status <- c(rexp(50, rate=1/30), rexp(50, rate=1/30))
status <- ifelse(time < status, 1, 0)
example_surv_1arm <- data.frame(status     = status,
                                time       = time,
                                historical = c(rep(1,50),rep(0,50)),
                                treatment  = 1)

fit01 <- bdpsurvival(Surv(time, status) ~ historical + treatment,
                     data = example_surv_1arm,
                     surv_time = 5, two_side=FALSE)
fit02 <- bdpsurvival(Surv(time, status) ~ historical + treatment,
                     data = example_surv_1arm,
                     surv_time = 5)

## ---- echo=FALSE---------------------------------------------------------
plot(fit02, type="discount")

## ---- echo=FALSE---------------------------------------------------------
plot(fit01, type="discount")

## ------------------------------------------------------------------------
p1 <- plot(fit01, type="discount", print=FALSE)
p1 + ggtitle("Discount Function Plot :-)")

## ------------------------------------------------------------------------
set.seed(42)
# Simulate survival times
time_current    <- rexp(50, rate=1/10)
time_historical <- rexp(50, rate=1/15)

# Combine simulated data into a data frame
data1 <- data.frame(status     = 1,
                    time       = c(time_current, time_historical),
                    historical = c(rep(0,50),rep(1,50)),
                    treatment  = 1)

## ------------------------------------------------------------------------
set.seed(42)
fit1 <- bdpsurvival(Surv(time, status) ~ historical + treatment,
                    data = data1,
                    surv_time = 5)

print(fit1)

## ---- include=FALSE------------------------------------------------------
survival_time_posterior_flat1 <- ppexp(5,
                                       fit1$posterior_treatment$posterior_hazard,
                                       cuts=c(0,fit1$args1$breaks))
surv_augmented1 <- round(1-median(survival_time_posterior_flat1), 4)
CI95_augmented1 <- round(1-quantile(survival_time_posterior_flat1, prob=c(0.975, 0.025)), 4)

## ------------------------------------------------------------------------
summary(fit1)

## ------------------------------------------------------------------------
set.seed(42)
fit1a <- bdpsurvival(Surv(time, status) ~ historical + treatment,
                     data = data1,
                     surv_time = 5,
                     alpha_max = 1,
                     fix_alpha = TRUE)

print(fit1a)

## ------------------------------------------------------------------------
survival_time_posterior_flat <- ppexp(5,
                                      fit1a$posterior_treatment$posterior_hazard,
                                      cuts=c(0,fit1a$args1$breaks))
surv_augmented <- 1-median(survival_time_posterior_flat)
CI95_augmented <- 1-quantile(survival_time_posterior_flat, prob=c(0.975, 0.025))

## ------------------------------------------------------------------------
plot(fit1a)

## ------------------------------------------------------------------------
set.seed(42)
# Simulate survival times for treatment data
time_current_trt    <- rexp(50, rate=1/10)
time_historical_trt <- rexp(50, rate=1/15)

# Simulate survival times for control data
time_current_cntrl    <- rexp(50, rate=1/20)
time_historical_cntrl <- rexp(50, rate=1/20)


# Combine simulated data into a data frame
data2 <- data.frame(status     = 1,
                    time       = c(time_current_trt,   time_historical_trt,
                                   time_current_cntrl, time_historical_cntrl),
                    historical = c(rep(0,50),rep(1,50), rep(0,50),rep(1,50)),
                    treatment  = c(rep(1,100), rep(0,100)))

## ------------------------------------------------------------------------
set.seed(42)
fit2 <- bdpsurvival(Surv(time, status) ~ historical + treatment,
                    data = data2)

print(fit2)

