context("bayesDP")

test_that("Normal", {
  expect_true(class(try(bdpnormal(20,5,500,30,10,250), silent = TRUE)) != "try-error",
               "1arm Normal Test Error")
  expect_true(class(try(bdpnormal(20,5,500,30,10,250,15,10,400,50,15,1000), silent = TRUE)) != "try-error",
               "2arm Normal Test/Control Error")
  options(warn=2)
  expect_true(class(try(bdpnormal(20,5,500,30,10,250), silent = TRUE)) != "try-error",
                 "1arm Normal Test Warning")
  expect_true(class(try(bdpnormal(9,9,9,9,9,9,9,9,9,9,9,9), silent = TRUE)) != "try-error",
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
  expect_true(class(try(bdpbinomial(30,200,20,700), silent = TRUE)) != "try-error",
              "1arm Binomial Test Error")
  expect_true(class(try(bdpbinomial(9,9,9,9,9,9,9,9,9,9,9,9), silent = TRUE)) != "try-error",
              "2arm Binomial Test/Control Error")
  options(warn=2)
  expect_true(class(try(bdpbinomial(30,200,20,700), silent = TRUE)) != "try-error",
                 "1arm Binomial Test Warning")
  expect_true(class(try(bdpbinomial(9,9,9,9,9,9,9,9,9,9,9,9), silent = TRUE)) != "try-error",
                 "2arm Binomial Test/Control Warning")
  options(warn=0)
})

test_that("Binomial Plot", {
  expect_true(class(try(plot(bdpbinomial(30,200,20,700)), silent = TRUE)) != "try-error",
              "1arm Binomial Test Plot Error")
  expect_true(class(try(plot(bdpbinomial(9,9,9,9,9,9,9,9)), silent = TRUE)) != "try-error",
               "2arm Binomial Test/Control Plot Error")
  options(warn=2)
  expect_true(class(try(plot(bdpbinomial(30,200,20,700)), silent = TRUE)) != "try-error",
                 "1arm Binomial Test Plot Warning")
  expect_true(class(try(plot(bdpbinomial(9,9,9,9,9,9,9,9)), silent = TRUE)) != "try-error",
                 "2arm Binomial Test/Control Plot Warning")
  options(warn=0)
})

test_that("Binomial Print", {
  expect_true(class(try(print(bdpbinomial(30,200,20,700)), silent = TRUE)) != "try-error",
              "1arm Binomial Test Print Error")
  expect_true(class(try(print(bdpbinomial(9,9,9,9,9,9,9,9)), silent = TRUE)) != "try-error",
               "2arm Binomial Test/Control Print Error")
  options(warn=2)
  expect_true(class(try(print(bdpbinomial(30,200,20,700)), silent = TRUE)) != "try-error",
                 "1arm Binomial Test Print Warning")
  expect_true(class(try(print(bdpbinomial(9,9,9,9,9,9,9,9)), silent = TRUE)) != "try-error",
                 "2arm Binomial Test/Control Print Warning")
  options(warn=0)
})

test_that("Binomial Summary", {
  expect_true(class(try(summary(bdpbinomial(30,200,20,700)), silent = TRUE)) != "try-error",
              "1arm Binomial Test Summary Error")
  expect_true(class(try(summary(bdpbinomial(9,9,9,9,9,9,9,9)), silent = TRUE)) != "try-error",
              "2arm Binomial Test/Control Summary Error")
  options(warn=2)
  expect_true(class(try(summary(bdpbinomial(30,200,20,700)), silent = TRUE)) != "try-error",
              "1arm Binomial Test Summary Warning")
  expect_true(class(try(summary(bdpbinomial(9,9,9,9,9,9,9,9)), silent = TRUE)) != "try-error",
              "2arm Binomial Test/Control Summary Warning")
  options(warn=0)
})
