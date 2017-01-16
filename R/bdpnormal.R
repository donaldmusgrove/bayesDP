#' bdpnormal
#'
#' bdpnormal
#'
#' @title bdpnormal: bdpnormal
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
#' @rdname bdpnormal
#' @export bdpnormal


setGeneric("bdpnormal",
function(type = c("1arm", "2arm"), ...){
  standardGeneric("bdpnormal")
})

setMethod("bdpnormal",
          signature(type = "character"),
          function(type = c("1arm", "2arm"), ...){
            if (type == "1arm"){
              return(bdpnormal1arm(...))
            }
            if(type == "2arm"){
              return(bdpnormal2arm(...))
            }
          })
