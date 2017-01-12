#' BayesianLossFunctionNormal
#'
#' BayesianLossFunctionNormal
#'
#' @title BayesianLossFunctionNormal: BayesianLossFunctionNormal
#' @param mu numeric
#' @param sigma2 numeric
#' @param N numeric
#' @param mu0 numeric
#' @param sigma02 numeric
#' @param type character
#' @param N0 numeric
#' @param alpha_max numeric
#' @param weibull_scale numeric
#' @param weibull_shape numeric
#' @param number_mcmc numeric
#' @param H0 numeric
#' @param two_side numeric
#' @param inequality character
#'
#' @examples
#'
#' @rdname BayesianLossFunctionNormal
#' @export BayesianLossFunctionNormal


setGeneric("BayesianLossFunctionNormal",
function(type = c("OPC", "RTC"), ...){
  standardGeneric("BayesianLossFunctionNormal")
})

setMethod("BayesianLossFunctionNormal",
          signature(type = "character"),
          function(type = c("OPC", "RTC"), ...){
            if (type == "OPC"){
              BayesianLossFunctionNormalOPC(...)
            }
            if(type == "RTC"){
              BayesianLossFunctionNormalRTC(...)
            }
          })
