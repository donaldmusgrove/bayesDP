################################################################################
# This code is used for estimating posterior samples from a Gaussian outcome   #
# where an informative prior is used. The prior weight is determined using a   #
# loss function. In addition this code simulate many trials in order to get    #
# trial characteristics you must specify the parameters of the loss function   #
# as well as the maximum strength for the prior. This code assumes a           #
# non-adaptive trial. This code is modeled after the methodologies developed   #
# by the MDIC working group:                                                   #
# 'Informing clinical trials using bench & simulations'                        #
# Developer: Tarek Haddad                                                      #
# Tarek.D.Haddad@Medtronic.com                                                 #
# Last modified: 9/10/2016                                                     #
################################################################################


################################################################################
### Results
################################################################################

## declare the display generic
setGeneric("rtcnormal", function(f, ...)
  standardGeneric("rtcnormal")
)

setMethod("rtcnormal",
          signature(f = "ANY"),
          function(f){
            message("Wrong object")
          })

setMethod("rtcnormal",
          signature(f = "missing"),
          function(f){
            message("Missing object")
          })

setMethod("rtcnormal",
          signature(f = "numeric"),
          function(
            loss_function_inputs = list(test_group_info = list(N0_max = 10,
                                                             weibull_scale = 0.1,
                                                             weibull_shape = 2),
                                        control_group_info = list(N0_max = 10,
                                                                  weibull_scale = 0.1,
                                                                  weibull_shape = 2)),

            test_group_info = list(test_group_info = list(mu = 10,
                                                          sigma = 0.1,
                                                          N = 2),
                                   control_group_info = list(mu = 10,
                                                             sigma = 0.1,
                                                             N = 2)),

            control_group_info = list(test_group_info = list(mu = 10,
                                                             sigma = 0.1,
                                                             N = 2),
                                      control_group_info = list(mu = 10,
                                                                sigma = 0.1,
                                                                N = 2)),
            two_sided_loss_function = False,
            number_mcmc = 1000
          ){



          })
