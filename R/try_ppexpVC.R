#' ppexpVC
#'
#' ppexpVC
#'
#' @title ppexpVC: ppexpVC
#' @rdname ppexpVC
#' @export ppexpVC
#' @importFrom Rcpp evalCpp

try_ppexpVC <- function() {
  .Call('bayesDP_ppexpVC', PACKAGE = 'bayesDP')
}

