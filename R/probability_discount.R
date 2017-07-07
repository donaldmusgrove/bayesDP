
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @title Bayesian Discount Prior: Comparison Between Current and Historical Data
#' @description \code{probability_discount} can be used to estimate the posterior
#'   probability of the comparison between historical and current data in the
#'   context of a one- or two-arm clinical trial with normal (mean) data.
#'   \code{probability_discount} is not used internally but is
#'   given for educational purposes.
#' @param mu scalar. Mean of the current data.
#' @param sigma scalar. Standard deviation of the current data.
#' @param N scalar. Number of observations of the current data.
#' @param mu0 scalar. Mean of the historical data.
#' @param sigma0 scalar. Standard deviation of the historical data.
#' @param N0 scalar. Number of observations of the historical data.
#' @param number_mcmc scalar. Number of Monte Carlo simulations. Default is 10000.
#' @param two_side logical. Indicator of two-sided test for the discount
#'   function. Default value is TRUE.
#' @details
#'   This function is not used internally but is given for educational purposes.
#'   Given the inputs,  the output is the posterior probability of the comparison
#'   between current and historical data in the context of a one- or two-arm clinical
#'   trial with normal (mean) data.
#'
#' @return \code{probability_discount} returns an object of class "probability_discount".
#'
#' An object of class \code{probability_discount} contains the following:
#' \describe{
#'  \item{\code{p_hat}}{
#'    scalar. The posterior probability of the comparison historical data weight.}
#' }
#'
#' @references
#' Haddad, T., Himes, A., Thompson, L., Irony, T., Nair, R. MDIC Computer
#'   Modeling and Simulation working group.(2017) Incorporation of stochastic
#'   engineering models as prior information in Bayesian medical device trials.
#'   \emph{Journal of Biopharmaceutical Statistics}, 1-15.
#'
#' @examples
#' probability_discount(mu  = 0,   sigma = 1, N  = 100,
#'                      mu0 = 0.1, sigma0 = 1, N0 = 100)
#'
#' @rdname probability_discount
#' @import methods
#' @importFrom stats sd density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov pchisq
#' @aliases probability_discount,ANY-method
#' @export probability_discount
probability_discount <- setClass("probability_discount")

setGeneric("probability_discount",
           function(mu          = NULL,
                    sigma       = NULL,
                    N           = NULL,
                    mu0         = NULL,
                    sigma0      = NULL,
                    N0          = NULL,
                    number_mcmc = 10000,
                    two_side    = TRUE){
             standardGeneric("probability_discount")
           })

setMethod("probability_discount",
          signature(),
          function(mu          = NULL,
                   sigma       = NULL,
                   N           = NULL,
                   mu0         = NULL,
                   sigma0      = NULL,
                   N0          = NULL,
                   number_mcmc = 10000,
                   two_side    = TRUE){

  ### Preposterior of current mu using flat prior
  sigma2_post_flat <- 1/rgamma(number_mcmc, (N - 1)/2, ((N - 1) * sigma^2)/2)
  mu_post_flat     <- rnorm(number_mcmc, mu, (sigma2_post_flat/((N-1)+1))^0.5)

  ### Posterior of historical data parameters using flat prior
  sigma2_post_flat0 <- 1/rgamma(number_mcmc, (N0-1)/2, ((N0-1)*sigma0^2)/2)
  mu_post_flat0     <- rnorm(number_mcmc, mu0, (sigma2_post_flat0/((N0-1)+1))^0.5)

  ### Test of model vs real
  p_hat <- mean(mu_post_flat < mu_post_flat0)

  ### Check if two-sided and transform p_hat accordingly
  if(two_side){
    p_hat <- ifelse(p_hat > 0.5, 1 - p_hat, p_hat)
  }

  return(p_hat)
})
