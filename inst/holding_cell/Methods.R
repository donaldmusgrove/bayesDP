#' plot
#'
#' plot
#'
#' @title plot: plot
#' @param x BayesianLossFunctionNormalOPC
#'
#' @examples
#'
#' @rdname plot
#' @export plot
setMethod("plot", signature(x = "BayesianLossFunctionNormalOPC"), function(x){
  op <- par(ask=TRUE)
  plot(x@post_typeplot1)
  plot(x@densityplot1)
  plot(x@lossfun_plot1)
  plot(x@lossfun_plot2)
  par(op)
})

#' plot
#'
#' plot
#'
#' @title print: print
#' @param x BayesianLossFunctionNormalOPC
#'
#' @examples
#'
#' @rdname print
#' @export print
setMethod("print", signature(x = "BayesianLossFunctionNormalOPC"), function(x){
  print(cat(x@hypothesis1))
  print(x@prior_for_test_group1)
})

