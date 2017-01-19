x <- airquality[, -1] # x is a regression design matrix
y <- airquality[,  1] # y is the corresponding response

stopifnot(complete.cases(y) != is.na(y))
ok <- complete.cases(x, y)
sum(!ok) # how many are not "ok" ?
x <- x[ok,]
y <- y[ok]

sapply(current,mean)
sapply(historical,sd)


#' data_helper
#'
#' data_helper
#'
#' @title data_helper: data_helper
#' @param mu_t numeric
#' @param sigma2_t numeric
#' @param N_t numeric
#' @param mu_c numeric
#' @param sigma2_c numeric
#' @param N_c numeric
#' @param mu0_t numeric
#' @param sigma02_t numeric
#' @param N0_t numeric
#' @param mu0_c numeric
#' @param sigma02_c numeric
#' @param N0_c numeric
#' @param type character
#' @param alpha_max numeric
#' @param weibull_scale numeric
#' @param weibull_shape numeric
#' @param number_mcmc numeric
#' @param two_side character
#' @param inequality character
#' @param delta numeric
#'
#'
#' @examples
#'
#' @rdname data_helper
#' @export data_helper


setGeneric("data_helper",
           function(current,historical){
             standardGeneric("data_helper")
           })

setMethod("data_helper",
          signature(historical = "missing"),
          function(type = c("1arm", "2arm"), ...){
            if (type == "1arm"){
              return(data_helper1arm(...))
            }
            if(type == "2arm"){
              return(data_helper2arm(...))
            }
          })

setMethod("data_helper",
          signature(historical = "missing"),
          function(type = c("1arm", "2arm"), ...){
            if (type == "1arm"){
              return(data_helper1arm(...))
            }
            if(type == "2arm"){
              return(data_helper2arm(...))
            }
          })


