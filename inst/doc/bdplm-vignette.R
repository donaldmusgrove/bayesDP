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
  
fit01 <- bdpbinomial(y_t=10, N_t=500, y0_t=25, N0_t=250)
fit02 <- bdpbinomial(y_t=10, N_t=500, y0_t=25, N0_t=250)

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
### Simulate  data
# Sample sizes
n_t  <- 25  #current treatment sample size
n_c  <- 25  #current control sample size
n_t0 <- 50  #historical treatment sample size
n_c0 <- 50  #historical treatment sample size

# Treatment and historical indicators
treatment  <- c(rep(1, n_t+n_t0), rep(0, n_c+n_c0))
historical <- c(rep(0, n_t), rep(1,n_t0), rep(0, n_c), rep(1,n_c0))

# Covariate effect
x    <- rnorm(n_t+n_c+n_t0+n_c0, 34, 5) 

# Outcome
Y  <- treatment + 0*historical + x*3.5 + rnorm(n_t+n_c+n_t0+n_c0,0,0.1)

# Place data in a single dataframe
df <- data.frame(Y=Y, treatment=treatment, historical=historical, x=x)

# Create current and historical dataframes
df_  <- subset(df, historical==0)
df_0 <- subset(df, historical==1)


# Fit the model with default inputs
fit <- bdplm(formula=Y ~ treatment+x,
             data=df_, data0=df_0)

# View estimates:
fit$estimates$coef

# View alpha discount weight parameters:
fit$alpha_discount


