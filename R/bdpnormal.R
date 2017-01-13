#' bdpnormal
#'
#' bdpnormal
#'
#' @title bdpnormal: bdpnormal
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
#' @rdname bdpnormal
#' @export bdpnormal


setGeneric("bdpnormal",
function(type = c("OPC", "RCT"), ...){
  standardGeneric("bdpnormal")
})

setMethod("bdpnormal",
          signature(type = "character"),
          function(type = c("OPC", "RCT"), ...){
            if (type == "OPC"){
              return(bdpnormal1arm(...))
            }
            if(type == "RCT"){
              return(bdpnormal2arm(...))
            }
          })
