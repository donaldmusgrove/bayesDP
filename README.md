![alt text](https://raw.githubusercontent.com/balcomes/bayesDP/master/bayesDP-logo.png "bayesDP Logo")

# bayesDP
### Tools for the Bayesian Discount Prior Function

https://cran.r-project.org/package=bayesDP

[![Travis-CI Build Status](https://travis-ci.org/balcomes/bayesDP.svg?branch=master)](https://travis-ci.org/balcomes/bayesDP)
[![Issue Count](https://codeclimate.com/github/balcomes/bayesDP/badges/issue_count.svg)](https://codeclimate.com/github/balcomes/bayesDP)
[![Downloads](http://cranlogs.r-pkg.org/badges/bayesDP?color=brightgreen)](http://www.r-pkg.org/pkg/bayesDP)
[![Project Status: Active ? The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

[![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/balcomes/bayesDP/issues)
[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%203%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![CRAN](http://www.r-pkg.org/badges/version/bayesDP)](https://cran.r-project.org/package=bayesDP)

### Description

Functions for data augmentation using the Bayesian discount prior function for 1 arm and 2 arm clinical trials.

### CRAN Installation

Install release version from CRAN:

```R
install.packages("bayesDP")
```

### GitHub Installation

Install development version from GitHub:

```R
devtools::install_github("balcomes/bayesDP")
```

### Documentation 

See manuals and vignettes within package.

### Examples

See manuals and vignettes within package.

{::nomarkdown}

<style type="text/css">
body, td {
   font-family: sans-serif;
   background-color: white;
   font-size: 13px;
}

body {
  max-width: 800px;
  margin: auto;
  padding: 1em;
  line-height: 20px;
}

tt, code, pre {
   font-family: 'DejaVu Sans Mono', 'Droid Sans Mono', 'Lucida Console', Consolas, Monaco, monospace;
}

h1 {
   font-size:2.2em;
}

h2 {
   font-size:1.8em;
}

h3 {
   font-size:1.4em;
}

h4 {
   font-size:1.0em;
}

h5 {
   font-size:0.9em;
}

h6 {
   font-size:0.8em;
}

a:visited {
   color: rgb(50%, 0%, 50%);
}

pre, img {
  max-width: 100%;
}
pre {
  overflow-x: auto;
}
pre code {
   display: block; padding: 0.5em;
}

code {
  font-size: 92%;
  border: 1px solid #ccc;
}

code[class] {
  background-color: #F8F8F8;
}

table, td, th {
  border: none;
}

blockquote {
   color:#666666;
   margin:0;
   padding-left: 1em;
   border-left: 0.5em #EEE solid;
}

hr {
   height: 0px;
   border-bottom: none;
   border-top-width: thin;
   border-top-style: dotted;
   border-top-color: #999999;
}

@media print {
   * {
      background: transparent !important;
      color: black !important;
      filter:none !important;
      -ms-filter: none !important;
   }

   body {
      font-size:12pt;
      max-width:100%;
   }

   a, a:visited {
      text-decoration: underline;
   }

   hr {
      visibility: hidden;
      page-break-before: always;
   }

   pre, blockquote {
      padding-right: 1em;
      page-break-inside: avoid;
   }

   tr, img {
      page-break-inside: avoid;
   }

   img {
      max-width: 100% !important;
   }

   @page :left {
      margin: 15mm 20mm 15mm 10mm;
   }

   @page :right {
      margin: 15mm 10mm 15mm 20mm;
   }

   p, h2, h3 {
      orphans: 3; widows: 3;
   }

   h2, h3 {
      page-break-after: avoid;
   }
}
</style>



</head>

<body>
<p>\name{bdpbinomial}
\alias{bdpbinomial}
\alias{bdpbinomial,ANY-method}
\title{Bayesian Discount Prior: Binomial counts}
\usage{
bdpbinomial(y_t = NULL, N_t = NULL, y0_t = NULL, N0_t = NULL,
  y_c = NULL, N_c = NULL, y0_c = NULL, N0_c = NULL, alpha_max = 1,
  fix_alpha = FALSE, a0 = 1, b0 = 1, number_mcmc = 10000,
  weibull_scale = 0.135, weibull_shape = 3, two_side = TRUE)
}
\arguments{
\item{y_t}{scalar. Number of events for the current treatment group.}</p>

<p>\item{N_t}{scalar. Sample size of the current treatment group.}</p>

<p>\item{y0_t}{scalar. Number of events for the historical treatment group.}</p>

<p>\item{N0_t}{scalar. Sample size of the historical treatment group.}</p>

<p>\item{y_c}{scalar. Number of events for the current control group.}</p>

<p>\item{N_c}{scalar. Sample size of the current control group.}</p>

<p>\item{y0_c}{scalar. Number of events for the historical control group.}</p>

<p>\item{N0_c}{scalar. Sample size of the historical control group.}</p>

<p>\item{alpha_max}{scalar. Maximum weight the discount function can apply.
Default is 1. For a two-arm trial, users may specify a vector of two values
where the first value is used to weight the historical treatment group and
the second value is used to weight the historical control group.}</p>

<p>\item{fix_alpha}{logical. Fix alpha at alpha_max? Default value is FALSE.}</p>

<p>\item{a0}{scalar. Prior value for the beta rate. Default is 1.}</p>

<p>\item{b0}{scalar. Prior value for the beta rate. Default is 1.}</p>

<p>\item{number_mcmc}{scalar. Number of Markov Chain Monte Carlo (MCMC)
simulations. Default is 10000.}</p>

<p>\item{weibull_scale}{scalar. Scale parameter of the Weibull discount function
used to compute alpha, the weight parameter of the historical data. Default
value is 0.135. For a two-arm trial, users may specify a vector of two values
where the first value is used to estimate the weight of the historical
treatment group and the second value is used to estimate the weight of the
historical control group.}</p>

<p>\item{weibull_shape}{scalar. Shape parameter of the Weibull discount function
used to compute alpha, the weight parameter of the historical data. Default
value is 3. For a two-arm trial, users may specify a vector of two values
where the first value is used to estimate the weight of the historical
treatment group and the second value is used to estimate the weight of the
historical control group.}</p>

<p>\item{two_side}{logical. Indicator of two-sided test for the discount
function. Default value is TRUE.}
}
\value{
\code{bdpbinomial} returns an object of class &ldquo;bdpbinomial&rdquo;.
The functions \code{summary} and \code{print} are used to obtain and
print a summary of the results, including user inputs. The \code{plot}
function displays visual outputs as well.</p>

<p>An object of class \code{bdpbinomial} is a list containing at least
the following components:</p>

<p>\describe{
 \item{\code{posterior_treatment}}{
   list. Entries contain values related to the treatment group:}
   \itemize{
     \item{\code{alpha_discount}}{
       numeric. Alpha value, the weighting parameter of the historical data.}
     \item{\code{p_hat}}{
       numeric. The posterior probability of the stochastic comparison
       between the current and historical data.}
     \item{\code{posterior}}{
       vector. The posterior of the treatment group, incorporating the
       weighted historical data.}
     \item{\code{posterior_flat}}{
       vector. The distribution of the current treatment group, i.e., no
       incorporation of the historical data.}
     \item{\code{prior}}{
       vector. The distribution of the historical treatment group.}
  }
 \item{\code{posterior_control}}{
   list. Similar entries as \code{posterior_treament}. Only present if
   control group is specified.}
 \item{\code{args1}}{
   list. Entries contain user inputs. In addition, the following elements
   are ouput:}
   \itemize{
     \item{\code{arm2}}{
       binary indicator. Used internally to indicate one-arm or two-arm
       analysis.}
     \item{\code{intent}}{
       character. Denotes current/historical status of treatment and
       control groups.}
  }
}
}
\description{
\code{bdpbinomial} is used for estimating posterior samples from a
  Binomial outcome where an informative prior is used. The prior weight
  is determined using a discount function. This code is modeled after
  the methodologies developed in Haddad et al. (2017).
}
\details{
\code{bdpbinomial} uses a two-stage approach for determining the
  strength of historical data in estimation of a binomial count mean outcome.
  In the first stage, a Weibull distribution function is used as a
  \emph{discount function} that defines the maximum strength of the
  historical data (via \code{weibull_shape}, \code{weibull_scale}, and
  \code{alpha_max}) and discounts based on disagreement with the current data.
  Disagreement between current and historical data is determined by stochastically
  comparing the respective posterior distributions under noninformative priors.
  With binomial data, the comparison is the proability (\code{p}) that the current
  count is less than the historical count. The comparison metric \code{p} is then
  input into the Weibull discount function and the final strength of the
  historical data is returned (alpha).</p>

<p>In the second stage, posterior estimation is performed where the discount
 function parameter, \code{alpha}, is used as a fixed value for all posterior
 estimation procedures.</p>

<p>To carry out a single arm (OPC) analysis, data for the current treatment
 (\code{y_t} and \code{N_t}) and historical treatment (\code{y0_t} and
 \code{N0_t}) must be input. The results are then based on the posterior
 distribution of the current data augmented by the historical data.</p>

<p>To carry our a two-arm (RCT) analysis, data for the current treatment and
 current control (\code{y_c} and \code{N_c}) must be input,
 as well as at least one of the historical treatment and historical control
 (\code{y0_c} and \code{N0_c}). The results
 are then based on the posterior distribution of the difference between
 current treatment and control, augmented by available historical data.
}
\examples{</p>

<h1>One-arm trial (OPC) example</h1>

<p>fit &lt;- bdpbinomial(y_t           = 10,
                   N_t           = 500,
                   y0_t          = 25,
                   N0_t          = 250)
summary(fit)
print(fit)
#plot(fit)</p>

<h1>Two-arm (RCT) example</h1>

<p>fit2 &lt;- bdpbinomial(y_t = 10,
                    N_t = 500,
                    y0_t = 25,
                    N0_t = 250,
                    y_c = 8,
                    N_c = 500,
                    y0_c = 20,
                    N0_c = 250)
summary(fit2)
print(fit2)
#plot(fit2)</p>

<p>}</p>

{:/}

### Authors

Shawn Balcome, Donnie Musgrove, Tarek Haddad
and Christopher Jackson (For the ppexp R code that was ported to C++.)

### License

GPL (>= 3)

:apple: :tangerine: :lemon: :cherries:  [:watermelon:](http://codeology.braintreepayments.com/balcomes/bayesdp#)  :strawberry:  :peach: :pear:  :green_apple:
