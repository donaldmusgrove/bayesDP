#' BayesianLossFunctionRTC
#'
#' BayesianLossFunctionRTC
#'
#' @title BayesianLossFunctionRTC: BayesianLossFunctionRTC
#' @param y numeric
#' @param N numeric
#' @param y0 numeric
#' @param N0 numeric
#' @param N0_max numeric
#' @param N0_test numeric
#' @param N0_control numeric
#' @param alpha0 numeric
#' @param beta0 numeric
#' @param number_mcmc numeric
#' @param weibull_shape numeric
#' @param weibull_scale numeric
#' @param two_side character
#'
#' @examples
#'
#' @rdname BayesianLossFunctionRTC
#' @export BayesianLossFunctionRTC


################################################################################
# This code is used for estimating posterior samples from a binary outcome     #
# where an informative prior is used. The prior weight is determined using a   #
# loss function. In addition this code simulate many trials in order to get    #
# trial characteristics you must specify the parameters of the loss function   #
# as well as the maximum strength for the prior. This code assumes a           #
# non-adaptive trial.                                                          #
# This code is modeled after the methodologies developed by the MDIC working   #
# group: "Informing clinical trials using bench & simulations"                 #
# Section 1: of the code defines functions needed                              #
# Section 2: of the code estimates a posterior and the loss function value     #
#            given inputs                                                      #
# Section 3: of the code simulates the trial many times to get trial           #
#            characteristics                                                   #
# Developer: Tarek Haddad                                                      #
# Tarek.D.Haddad@Medtronic.com                                                 #
# Last modified:1/26/2016                                                      #
################################################################################
#library(ggplot2)

setGeneric("BayesianLossFunctionRTC",
           function(y             = 1,          #Number of events observed  current data sets
                    N             = 400,        #Number of  current subjects
                    y0            = 10,         #Number of events observed  historical  data sets
                    N0            = 100, #Number of historical subjects
                    N0_max        = 200,        #Max effective sample size the prior can receive when loss function equals 1
                    N0_test       = 1000,
                    N0_control    = 100,
                    alpha0        = 1,          #Noninformative Initial priors
                    beta0         = 1,          #Noninformative Initial priors
                    number_mcmc   = 10000,      #Number of simulations to estimate posterior and loss function
                    weibull_scale = .05,        #Loss function parameter controlling the location of a weibull function
                    weibull_shape = 2,          #Loss function parameter controlling the location of a weibull function
                    two_side      = 0){
             standardGeneric("BayesianLossFunctionRTC")
           })

setMethod("BayesianLossFunctionRTC",
          signature(y = "numeric"),
          function(y             = 1,          #Number of events observed  current data sets
                   N             = 400,        #Number of  current subjects
                   y0            = 10,         #Number of events observed  historical  data sets
                   N0            = 100, #Number of historical subjects
                   N0_max        = 200,        #Max effective sample size the prior can receive when loss function equals 1
                   N0_test       = 1000,
                   N0_control    = 100,
                   alpha0        = 1,          #Noninformative Initial priors
                   beta0         = 1,          #Noninformative Initial priors
                   number_mcmc   = 10000,      #Number of simulations to estimate posterior and loss function
                   weibull_scale = .05,        #Loss function parameter controlling the location of a weibull function
                   weibull_shape = 2,          #Loss function parameter controlling the location of a weibull function
                   two_side      = 0){


### Section 1: Functions
############################################################
#This function produces the appropriate weight for your    #
#prior data assuming a binomial outcome                    #
#This function produces a scaler between 0 and 1           #
############################################################
Loss_function <- function(y,N,y0,N0,alpha0,beta0,number_mcmc,weibull_shape,
                          weibull_scale, two_side=0){

  ### Theta for using Flat prior
  alpha_post_flat  <- y+alpha0
  beta_post_flat   <- N-y+beta0
  theta_post_flate <- rbeta(number_mcmc, alpha_post_flat, beta_post_flat)

  ### Prior model###########
  alpha_prior <- y0 +alpha0
  beta_prior  <- N0-y0+beta0
  theta0      <- rbeta(number_mcmc,  alpha_prior,beta_prior)#  flate prior

  ### Test of model vs real
  p_test <- mean(theta_post_flate<theta0)   # larger is higher failure

  ### Number of effective sample size given shape and scale loss function
  if (two_side == 0) {
    alpha_loss <- pweibull(p_test, shape = weibull_shape, scale = weibull_scale)
  } else if (two_side == 1){
    p_test1    <- ifelse(p_test > 0.5, 1 - p_test, p_test)
    alpha_loss <- pweibull(p_test1, shape = weibull_shape, scale = weibull_scale)
  }

  return(list(alpha_loss       = alpha_loss,
              pvalue           = p_test,
              theta_post_flate = theta_post_flate,
              theta0           = theta0))
}

#############################################################
# Calculates posterior estimation for Binomial distribution #
# Given alpha_loss value and maximum strength of(N0_max)    #
# prior if alpha_loss = 1                                   #
#############################################################
theta_post_aug_bin <- function(y,N,y0,N0,N0_max,alpha_loss,alpha0,beta0,
                               number_mcmc){

  effective_N0 <- N0_max*alpha_loss

  if(N0==0){
    alpha_prior <- alpha0
    beta_prior  <- beta0
  }else{
    alpha_prior <- (y0/N0)*effective_N0+alpha0
    beta_prior  <- effective_N0-(y0/N0)*effective_N0+beta0
  }

  alpha_post_aug <- y+alpha_prior
  beta_post_aug  <- N-y+beta_prior
  theta_post_aug <- rbeta(number_mcmc,alpha_post_aug,beta_post_aug)
  return(theta_post_aug)
}

############################################################
# Combines the loss function and posterior                 #
# estimation into one function                             #
############################################################
Binomial_posterior <- function(y,N,y0,N0,N0_max,alpha0,beta0,number_mcmc,
                               weibull_shape,weibull_scale,two_side){

  alpha_loss         <- Loss_function(y=y,N=N,y0=y0,N0=N0,alpha0=alpha0,
                                      beta0=beta0,number_mcmc=number_mcmc,
                                      weibull_shape=weibull_shape,
                                      weibull_scale=weibull_scale,
                                      two_side=two_side)
  Binomial_posterior <- theta_post_aug_bin(y=y,N=N,y0=y0,N0=N0,N0_max=N0_max,
                                           alpha_loss=alpha_loss$alpha_loss,
                                           alpha0=alpha0,
                                           beta0=beta0,number_mcmc=number_mcmc)

  return(list(alpha_loss               = alpha_loss$alpha_loss,
              pvalue                   = alpha_loss$pvalue,
              Binomial_posterior       = Binomial_posterior,
              Binomial_posterior_flate = alpha_loss$theta_post_flate,
              Binomial_prior           = alpha_loss$theta0,
              weibull_scale            = weibull_scale,
              weibull_shape            = weibull_shape,
              y                        = y,
              N                        = N,
              y0                       = y0,
              N0                       = N0,
              N0_max                   = N0_max,
              N0_effective             = alpha_loss$alpha_loss*N0_max))
}


final <- function(posterior_control,posterior_test){
  den_post_control  <- density(posterior_control$Binomial_posterior,adjust = .5)
  den_flat_control  <- density(posterior_control$Binomial_posterior_flate,adjust = .5)
  den_prior_control <- density(posterior_control$Binomial_prior,adjust = .5)

  den_post_test     <- density(posterior_test$Binomial_posterior,adjust = .5)
  den_flat_test     <- density(posterior_test$Binomial_posterior_flate,adjust = .5)
  den_prior_test    <- density(posterior_test$Binomial_prior,adjust = .5)

  TestMinusControl_post <- posterior_test$Binomial_posterior-posterior_control$Binomial_posterior

  return(list(den_post_control      = den_post_control,
              den_flat_control      = den_flat_control,
              den_prior_control     = den_prior_control,
              den_post_test         = den_post_test,
              den_flat_test         = den_flat_test,
              den_prior_test        = den_prior_test,
              TestMinusControl_post = TestMinusControl_post))
}




#N0_test    <- 1000
#N0_control <- 100

### Function call & Input Variable name description
### - Input Variable name description
posterior_test <- Binomial_posterior(
  y,             #= 10,      #Number of events observed  current data sets
  N,             #= 400,     #Number of  current subjects
  y0,            #= 1,       #Number of events observed  historical  data sets
  N0,            #= N0_test, #Number of historical subjects
  N0_max,        #= 200,     #Max effective sample size the prior can receive when loss function equals 1
  alpha0,        #= 1,       #Noninformative Initial priors
  beta0,         #= 1,       #Noninformative Initial priors
  number_mcmc,   #= 10000,   #Number of simulations to estimate posterior and loss function
  weibull_scale, #= .2,      #Loss function parameter controlling the location of a weibull function
  weibull_shape, #= 2,       #Loss function parameter controlling the location of a weibull function
  two_side       #= 0
)

posterior_control <- Binomial_posterior(
  y,             #= 1,          #Number of events observed  current data sets
  N,             #= 400,        #Number of  current subjects
  y0,            #= 10,         #Number of events observed  historical  data sets
  N0,            #= N0_control, #Number of historical subjects
  N0_max,        #= 200,        #Max effective sample size the prior can receive when loss function equals 1
  alpha0,        #= 1,          #Noninformative Initial priors
  beta0,         #= 1,          #Noninformative Initial priors
  number_mcmc,   #= 10000,      #Number of simulations to estimate posterior and loss function
  weibull_scale, #= .05,        #Loss function parameter controlling the location of a weibull function
  weibull_shape, #= 2,          #Loss function parameter controlling the location of a weibull function
  two_side       #= 0
)

f <- final(posterior_control,posterior_test)



### Results
### Print

if(posterior_test$N0==0){
  print("No Prior Supplied")
}else{
  list("Effective sample size of prior"=posterior_test$N0_effective,
       "Bayesian p-value (new vs historical data)"=posterior_test$pvalue,
       'loss function value'=posterior_test$alpha_loss,
       'N0_max'=posterior_test$N0_max
  )
}

if(posterior_control$N0==0){
  print("No Prior Supplied")
}else{
  list("Effective sample size of prior"=posterior_control$N0_effective,
       "Bayesian p-value (new vs historical data)"=posterior_control$pvalue,
       'loss function value'=posterior_control$alpha_loss,
       'N0_max'=posterior_control$N0_max
  )

}


D1 <- data.frame(information_sources='Posterior',group="Control",y=f$den_post_control$y,x=f$den_post_control$x)
D2 <- data.frame(information_sources="Current data",group="Control",y=f$den_flat_control$y,x=f$den_flat_control$x)
D3 <- data.frame(information_sources="Prior",group="Control",y=f$den_prior_control$y,x=f$den_prior_control$x)
D4 <- data.frame(information_sources='Posterior',group="Test",y=f$den_post_test$y,x=f$den_post_test$x)
D5 <- data.frame(information_sources="Current data",group="Test",y=f$den_flat_test$y,x=f$den_flat_test$x)
D6 <- data.frame(information_sources="Prior",group="Test",y=f$den_prior_test$y,x=f$den_prior_test$x)

if(N0_test==0 & N0_control==0){
  D <- as.data.frame(rbind(D1,D2,D4,D5))
}
if(N0_test==0 & N0_control!=0){
  D <- as.data.frame(rbind(D1,D2,D3,D4,D5))
}
if(N0_test!=0 & N0_control==0){
  D <- as.data.frame(rbind(D1,D2,D4,D5,D6))
}
if(N0_test!=0 & N0_control!=0){
  D <- as.data.frame(rbind(D1,D2,D3,D4,D5,D6,D6))
}

D$information_sources <- factor(D$information_sources,levels = (c("Posterior","Current data","Prior")))

p <- ggplot(D,aes(x=x,y=y)) +
  geom_line(size=2,aes(colour=information_sources,lty=information_sources)) +
  theme_bw() +
  facet_wrap(~group, ncol=1) +
  ylab("Density (PDF)") +
  xlab("values")


p1 <- ggplot(subset(D,information_sources=="Posterior"), aes(x=x,y=y)) +
  geom_line(size=2,aes(colour=group)) +
  ylab("Density (PDF)") +
  xlab("values") +
  theme_bw()


pHist <- ggplot(data.frame(x = f$TestMinusControl_post),
                aes(x=x))  +
  geom_histogram(fill="dodgerblue", color="white") +
  xlab("Test-Control") +
  geom_vline(xintercept = 0, linetype = "longdash")

p_value <- seq(0,1,,100)
Loss_function_test <- pweibull(p_value,
                               shape=posterior_test$weibull_shape,
                               scale=posterior_test$weibull_scale)*posterior_test$N0_max
Loss_function_control <- pweibull(p_value,
                                  shape=posterior_control$weibull_shape,
                                  scale=posterior_control$weibull_scale)*posterior_control$N0_max


D1 <- data.frame(group="test",y=Loss_function_test,x=p_value)
D2 <- data.frame(group=c("test"),pvalue=c(posterior_test$pvalue))
D3 <- data.frame(group=c("test"),pvalue=c(posterior_test$N0_effective))
D4 <- data.frame(group="control",y=Loss_function_control,x=p_value)
D5 <- data.frame(group=c("control"),pvalue=c(posterior_control$pvalue))
D6 <- data.frame(group=c("control"),pvalue=c(posterior_control$N0_effective))


p3 <- ggplot()
if(N0_test!=0){p3=p3+geom_line(data=D1,aes(y=y,x=x,colour=group),size=1)+
  geom_vline(data=D2, aes(xintercept=pvalue,colour=group),lty=2)+
  geom_hline(data=D3, aes(yintercept=pvalue,colour=group),lty=2)
}
if(N0_control!=0){p3=p3+geom_line(data=D4,aes(y=y,x=x,colour=group),size=1)+
  geom_vline(data=D5, aes(xintercept=pvalue,colour=group),lty=2)+
  geom_hline(data=D6, aes(yintercept=pvalue,colour=group),lty=2)
}
p3 <- p3+facet_wrap(~group, ncol=1) +
  theme_bw() +
  ylab("Effective sample size for historical data") +
  xlab("Bayesian p-value (new vs historical data)")


cat(paste('"We can define W as the difference between the event rates for the test group versus control group, i.e. W=test rate - control rate',
          '\n',
          'Null Hypothesis (H_0): W≥0.0%',
          '\n',
          "Alternative Hypothesis (H_a): W<0.0%",
          '\n',
          '\n',
          "P(W<0.0%|data)=",mean(f$TestMinusControl_post<0),
          '\n',
          "We can accept H_a with a Probability of",mean(f$TestMinusControl_post<0)
))

me <- list(cat(paste('"We can define W as the difference between the event rates for the test group versus control group, i.e. W=test rate - control rate',
                     '\n',
                     'Null Hypothesis (H_0): W≥0.0%',
                     '\n',
                     "Alternative Hypothesis (H_a): W<0.0%",
                     '\n',
                     '\n',
                     "P(W<0.0%|data)=",mean(f$TestMinusControl_post<0),
                     '\n',
                     "We can accept H_a with a Probability of",mean(f$TestMinusControl_post<0))),
           p,
           p1,
           p3,
           pHist)

## Set the name for the class
class(me) <- append(class(me),"BayesianLossFunctionRCT")
return(me)

})

