context("bayesDP")


################################################################################
# bdpnormal
################################################################################

# One-arm trial (OPC) example
fit <- bdpnormal(mu_t = 30, sigma_t = 10, N_t = 250,
                 mu0_t = 50, sigma0_t = 5, N0_t = 250)
summary(fit)
plot(fit)

# Two-arm (RCT) example
fit2 <- bdpnormal(mu_t = 30, sigma_t = 10, N_t = 250,
                  mu0_t = 50, sigma0_t = 5, N0_t = 250,
                  mu_c = 25, sigma_c = 10, N_c = 250,
                  mu0_c = 50, sigma0_c = 5, N0_c = 250)
summary(fit2)
plot(fit2)


################################################################################
# bdpbinomial
################################################################################


# One-arm trial (OPC) example
fit <- bdpbinomial(y_t           = 10,
                   N_t           = 500,
                   y0_t          = 25,
                   N0_t          = 250)
summary(fit)
print(fit)
plot(fit)

# Two-arm (RCT) example
fit2 <- bdpbinomial(y_t = 10,
                    N_t = 500,
                    y0_t = 25,
                    N0_t = 250,
                    y_c = 8,
                    N_c = 500,
                    y0_c = 20,
                    N0_c = 250)
summary(fit2)
print(fit2)
plot(fit2)


################################################################################
# bdpsurvival
################################################################################


# One-arm trial (OPC) example - survival probability at 5 years
# Simulate survival data for a single arm (OPC) trial
time   <- c(rexp(50, rate=1/20), rexp(50, rate=1/10))
status <- c(rexp(50, rate=1/30), rexp(50, rate=1/30))
status <- ifelse(time < status, 1, 0)

# Collect data into a data frame
example_surv_1arm <- data.frame(status     = status,
                                time       = time,
                                historical = c(rep(1,50),rep(0,50)),
                                treatment  = 1)

fit1 <- bdpsurvival(Surv(time, status) ~ historical + treatment,
                    data = example_surv_1arm,
                    surv_time = 5)

print(fit1)

plot(fit1)


# Two-arm trial (OPC) example
# Simulate survival data for a two-arm trial
time   <- c(rexp(50, rate=1/20), # Current treatment
            rexp(50, rate=1/10), # Current control
            rexp(50, rate=1/30), # Historical treatment
            rexp(50, rate=1/5))  # Historical control
status <- rexp(200, rate=1/40)
status <- ifelse(time < status, 1, 0)

# Collect data into a data frame
example_surv_2arm <- data.frame(status     = status,
                                time       = time,
                                historical = c(rep(0,100),rep(1,100)),
                                treatment  = c(rep(1,50),rep(0,50),rep(1,50),rep(0,50)))

fit2 <- bdpsurvival(Surv(time, status) ~ historical + treatment,
                    data = example_surv_2arm)

summary(fit2)

### Fix alpha at 1
fit2_1 <- bdpsurvival(Surv(time, status) ~ historical + treatment,
                      data = example_surv_2arm,
                      fix_alpha = TRUE)

summary(fit2_1)
