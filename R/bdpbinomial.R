#' bdpbinomial
#'
#' bdpbinomial
#'
#' @title bdpbinomial: bdpbinomial
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
