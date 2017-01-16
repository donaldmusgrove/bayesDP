#' bdpbinomial
#'
#' bdpbinomial
#'
#' @title bdpbinomial: bdpbinomial
#' @param y_t numeric
#' @param N_t numeric
#' @param y_c numeric
#' @param N_c numeric
#' @param y0_t numeric
#' @param N0_t numeric
#' @param y0_c numeric
#' @param N0_c numeric
#' @param type numeric
#' @param alpha_max numeric
#' @param a0 numeric
#' @param b0 numeric
#' @param number_mcmc numeric
#' @param weibull_shape numeric
#' @param weibull_scale numeric
#' @param two_side character
#' @param inequality character
#' @param delta character
#'
#' @examples
#'
#' @rdname bdpbinomial
#' @export bdpbinomial


setGeneric("bdpbinomial",
           function(type = c("1arm", "2arm"), ...){
             standardGeneric("bdpbinomial")
           })

setMethod("bdpbinomial",
          signature(type = "character"),
          function(type = c("1arm", "2arm"), ...){
            if (type == "1arm"){
              return(bdpbinomial1arm(...))
            }
            if(type == "2arm"){
              return(bdpbinomial2arm(...))
            }
          })
