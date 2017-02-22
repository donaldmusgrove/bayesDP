context("bayesDP")

test_that("Normal", {
  expect_true(class(try(bdpnormal(20,5,500,30,10,250), silent = TRUE)) != "try-error",
              "1arm Normal Test Error")
  expect_true(class(try(bdpnormal(20,5,500,30,10,250,15,10,400,50,15,1000), silent = TRUE)) != "try-error",
              "2arm Normal Test/Control Error")
  options(warn=2)
  expect_true(class(try(bdpnormal(20,5,500,30,10,250), silent = TRUE)) != "try-error",
              "1arm Normal Test Warning")
  expect_true(class(try(bdpnormal(20,5,500,30,10,250,15,10,400,50,15,1000), silent = TRUE)) != "try-error",
              "2arm Normal Test/Control Warning")
  options(warn=0)
})

test_that("Normal Plot", {
  expect_true(class(try(plot(bdpnormal(20,5,500,30,10,250)), silent = TRUE)) != "try-error",
              "1arm Normal Test Plot Error")
  expect_true(class(try(plot(bdpnormal(20,5,500,30,10,250,15,10,400,50,15,1000)), silent = TRUE)) != "try-error",
              "2arm Normal Test/Control Plot Error")
  options(warn=2)
  expect_true(class(try(plot(bdpnormal(20,5,500,30,10,250)), silent = TRUE)) != "try-error",
              "1arm Normal Test Plot Warning")
  expect_true(class(try(plot(bdpnormal(20,5,500,30,10,250,15,10,400,50,15,1000)), silent = TRUE)) != "try-error",
              "2arm Normal Test/Control Plot Warning")
  options(warn=0)
})

test_that("Normal Print", {
  expect_true(class(try(print(bdpnormal(20,5,500,30,10,250)), silent = TRUE)) != "try-error",
              "1arm Normal Test Print Error")
  expect_true(class(try(print(bdpnormal(20,5,500,30,10,250,15,10,400,50,15,1000)), silent = TRUE)) != "try-error",
              "2arm Normal Test/Control Print Error")
  options(warn=2)
  expect_true(class(try(print(bdpnormal(20,5,500,30,10,250)), silent = TRUE)) != "try-error",
              "1arm Normal Test Print Warning")
  expect_true(class(try(print(bdpnormal(20,5,500,30,10,250,15,10,400,50,15,1000)), silent = TRUE)) != "try-error",
              "2arm Normal Test/Control Print Warning")
  options(warn=0)
})

test_that("Normal Summary", {
  expect_true(class(try(summary(bdpnormal(20,5,500,30,10,250)), silent = TRUE)) != "try-error",
              "1arm Normal Test Summary Error")
  expect_true(class(try(summary(bdpnormal(20,5,500,30,10,250,15,10,400,50,15,1000)), silent = TRUE)) != "try-error",
              "2arm Normal Test/Control Summary Error")
  options(warn=2)
  expect_true(class(try(summary(bdpnormal(20,5,500,30,10,250)), silent = TRUE)) != "try-error",
              "1arm Normal Test Summary Warning")
  expect_true(class(try(summary(bdpnormal(20,5,500,30,10,250,15,10,400,50,15,1000)), silent = TRUE)) != "try-error",
              "2arm Normal Test/Control Summary Warning")
  options(warn=0)
})



test_that("Binomial", {
  expect_true(class(try(bdpbinomial(20,500,50,300), silent = TRUE)) != "try-error",
              "1arm Binomial Test Error")
  expect_true(class(try(bdpbinomial(20,500,50,300,30,400,60,1000), silent = TRUE)) != "try-error",
              "2arm Binomial Test/Control Error")
  options(warn=2)
  expect_true(class(try(bdpbinomial(20,500,50,300), silent = TRUE)) != "try-error",
              "1arm Binomial Test Warning")
  expect_true(class(try(bdpbinomial(20,500,50,300,30,400,60,1000), silent = TRUE)) != "try-error",
              "2arm Binomial Test/Control Warning")
  options(warn=0)
})

test_that("Binomial Plot", {
  expect_true(class(try(plot(bdpbinomial(20,500,50,300)), silent = TRUE)) != "try-error",
              "1arm Binomial Test Plot Error")
  expect_true(class(try(plot(bdpbinomial(20,500,50,300,30,400,60,1000)), silent = TRUE)) != "try-error",
              "2arm Binomial Test/Control Plot Error")
  options(warn=2)
  expect_true(class(try(plot(bdpbinomial(20,500,50,300)), silent = TRUE)) != "try-error",
              "1arm Binomial Test Plot Warning")
  expect_true(class(try(plot(bdpbinomial(20,500,50,300,30,400,60,1000)), silent = TRUE)) != "try-error",
              "2arm Binomial Test/Control Plot Warning")
  options(warn=0)
})

test_that("Binomial Print", {
  expect_true(class(try(print(bdpbinomial(20,500,50,300)), silent = TRUE)) != "try-error",
              "1arm Binomial Test Print Error")
  expect_true(class(try(print(bdpbinomial(20,500,50,300,30,400,60,1000)), silent = TRUE)) != "try-error",
              "2arm Binomial Test/Control Print Error")
  options(warn=2)
  expect_true(class(try(print(bdpbinomial(20,500,50,300)), silent = TRUE)) != "try-error",
              "1arm Binomial Test Print Warning")
  expect_true(class(try(print(bdpbinomial(20,500,50,300,30,400,60,1000)), silent = TRUE)) != "try-error",
              "2arm Binomial Test/Control Print Warning")
  options(warn=0)
})

test_that("Binomial Summary", {
  expect_true(class(try(summary(bdpbinomial(20,500,50,300)), silent = TRUE)) != "try-error",
              "1arm Binomial Test Summary Error")
  expect_true(class(try(summary(bdpbinomial(20,500,50,300,30,400,60,1000)), silent = TRUE)) != "try-error",
              "2arm Binomial Test/Control Summary Error")
  options(warn=2)
  expect_true(class(try(summary(bdpbinomial(20,500,50,300)), silent = TRUE)) != "try-error",
              "1arm Binomial Test Summary Warning")
  expect_true(class(try(summary(bdpbinomial(20,500,50,300,30,400,60,1000)), silent = TRUE)) != "try-error",
              "2arm Binomial Test/Control Summary Warning")
  options(warn=0)
})

test_that("Linear Regression", {
  data <- data.frame(y         = rnorm(100, 4, 0.1),
                     x         = c(rnorm(50,1,0.1), rnorm(50,3,0.1)),
                     treatment = c(rep(0,50),rep(1,50)))

  expect_true(class(try(bdpregression_linear(data,
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
                                             two_side      = 0), silent = TRUE)) != "try-error",
              "Linear Regression Error")
})
