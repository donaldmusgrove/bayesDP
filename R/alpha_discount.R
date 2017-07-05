
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @title INSERT TEXT HERE
#' @description \code{alpha_discount} INSERT TEXT HERE
#' @param p_hat INSERT TEXT HERE
#' @param alpha_max scalar. Maximum weight the discount function can apply.
#'   Default is 1. For a two-arm trial, users may specify a vector of two values
#'   where the first value is used to weight the historical treatment group and
#'   the second value is used to weight the historical control group.
#' @param weibull_shape scalar. Shape parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 3. For a two-arm trial, users may specify a vector of two values
#'   where the first value is used to estimate the weight of the historical
#'   treatment group and the second value is used to estimate the weight of the
#'   historical control group.
#' @param weibull_scale scalar. Scale parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 0.135. For a two-arm trial, users may specify a vector of two values
#'   where the first value is used to estimate the weight of the historical
#'   treatment group and the second value is used to estimate the weight of the
#'   historical control group.
#' @details \code{alpha_discount} INSERT TEXT HERE
#'
#'   For more details, see the \code{alpha_discount} vignette: \cr
#'   \code{vignette("alpha_discount-vignette", package="bayesDP")}
#'
#'
#' @return \code{alpha_discount} returns an object of class "alpha_discount".
#'
#' An object of class \code{alpha_discount} is a list containing at least
#' the following components:
#' \describe{
#'  \item{\code{posterior_treatment}}{
#'    list. Entries contain values related to the data:}
#'    \itemize{
#'     INSERT TEXT HERE
#'   }
#' }
#'
#' @seealso \code{\link[=summary,alpha_discount-method]{summary}},
#'   \code{\link[=print,alpha_discount-method]{print}},
#'   and \code{\link[=plot,alpha_discount-method]{plot}} for details of each of the
#'   supported methods.
#'
#' @references
#' Haddad, T., Himes, A., Thompson, L., Irony, T., Nair, R. MDIC Computer
#'   Modeling and Simulation working group.(2017) Incorporation of stochastic
#'   engineering models as prior information in Bayesian medical device trials.
#'   \emph{Journal of Biopharmaceutical Statistics}, 1-15.
#'
#' @examples
#' # INSERT EXAMPLE HERE
#'
#' @rdname alpha_discount
#' @import methods
#' @importFrom stats sd density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov pchisq
#' @export alpha_discount
alpha_discount <- setClass("alpha_discount")

setGeneric("alpha_discount",
           function(p_hat          = NULL,
                    alpha_max     = 1,
                    weibull_scale = 0.135,
                    weibull_shape = 3){
             standardGeneric("alpha_discount")
           })

setMethod("alpha_discount",
          signature(),
          function(p_hat          = NULL,
                   alpha_max     = 1,
                   weibull_scale = 0.135,
                   weibull_shape = 3){

  alpha_hat <- pweibull(p_hat, shape = weibull_shape, scale = weibull_scale) * alpha_max
  return(alpha_hat)
})
