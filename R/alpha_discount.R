
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @title Bayesian Discount Prior: Historical Data Weight (alpha)
#' @description \code{alpha_discount} can be used to estimate the weight
#'   applied to historical data in the context of a one- or two-arm
#'   clinical trial. \code{alpha_discount} is not used internally but is
#'   given for educational purposes.
#' @param p_hat scalar. The posterior probability of a stochastic comparison.
#'   This value can be the output of \code{posterior_probability} or a value
#'   between 0 and 1.
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
#' @details
#'   This function is not used internally but is given for educational purposes.
#'   Given inputs \code{p_hat}, \code{alpha_max}, \code{weibull_shape}, and
#'   \code{weibull_scale} the output is the weight that would be applied to
#'   historical data in the context of a one- or two-arm clinical trial.
#'
#' @return \code{alpha_discount} returns an object of class "alpha_discount".
#'
#' An object of class \code{alpha_discount} contains the following:
#' \describe{
#'  \item{\code{alpha_hat}}{
#'    scalar. The historical data weight.
#'   }
#'  }
#'
#' @references
#' Haddad, T., Himes, A., Thompson, L., Irony, T., Nair, R. MDIC Computer
#'   Modeling and Simulation working group.(2017) Incorporation of stochastic
#'   engineering models as prior information in Bayesian medical device trials.
#'   \emph{Journal of Biopharmaceutical Statistics}, 1-15.
#'
#' @examples
#' alpha_discount(0.5)
#'
#' alpha_discount(0)
#'
#' @rdname alpha_discount
#' @import methods
#' @importFrom stats sd density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov pchisq
#' @aliases alpha_discount,ANY-method
#' @export alpha_discount
alpha_discount <- setClass("alpha_discount")

setGeneric("alpha_discount",
           function(p_hat         = NULL,
                    alpha_max     = 1,
                    weibull_scale = 0.135,
                    weibull_shape = 3){
             standardGeneric("alpha_discount")
           })

setMethod("alpha_discount",
          signature(),
          function(p_hat         = NULL,
                   alpha_max     = 1,
                   weibull_scale = 0.135,
                   weibull_shape = 3){

  alpha_hat <- pweibull(p_hat, shape = weibull_shape, scale = weibull_scale) * alpha_max
  return(alpha_hat)
})
