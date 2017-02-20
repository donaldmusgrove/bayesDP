#' @title Compute cdf of the piecewise exponential distribution
#'
#' @Description Returns the cumulative distribution function computed at time \code{q},
#' hazard(s) \code{x}, and cut points \code{cuts}.
#'
#' @param q scalar. The time point at which the cdf is to be computed.
#' @param x vector of hazard rates.
#' @param cuts vector of the same length as x, giving the times at which the
#' rate changes, i.e., the interval cut points. The first element of cuts
#' should be 0, and cuts should be in increasing order.
#'
#' @details
#'
#' @examples
#' q    <- 12
#' x    <- c(0.25,0.3,0.35,0.4)
#' cuts <- c(0,6,12,18)
#' pp   <- ppexp(q,x,cuts)
#'
#' @useDynLib bayesDP
#'
#' @export
ppexp <- function(q, x, cuts){

  if(!is.matrix(x)){
    ppout <- ppexpV(q, x, cuts)
  } else if(is.matrix(x)){
    ppout <- ppexpM(q, x, cuts)
  } else{
    return("Error: input x is in the wrong format.")
  }

  return(ppout)
}
