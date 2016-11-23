########################################################################################
########################################################################################
#This code is used for estimating posterior samples from a gaussian outcome where an   #
#informative prior is used. The prior weight is determined using a loss function.      #
#In addition this code simulate many trials in order to get trial characteristics      #
#you must specify the parameters of the loss function as well as the maximum strength  #
#for the prior. This code assumes a non-adaptive trial                                 #
#This code is modeled after the methodologies developed by the MDIC working group:     #
#"Informing clinical trials using bench & simulations"                                 #
#Section 1:  of the code defines functions needed                                      #
#Section 2:  of the code estimates a posterior and the loss function value given inputs#
#Section 3:  of the code simulates the trial many times to get trial characteristics   #
#Developer: Tarek Haddad                                                               #
#Tarek.D.Haddad@Medtronic.com                                                          #
#Last modified:1/26/2016                                                               #


################################################################################
### Results     
################################################################################

results <- function(f,posterior_test,H0,two_side,inequality){

  D4 <- data.frame(information_sources = "Posterior", 
                   group               = "Test", 
                   y                   = f$den_post_test$y, 
                   x                   = f$den_post_test$x)
                     
  D5 <- data.frame(information_sources = "Current data", 
                   group               = "Test", 
                   y                   = f$den_flat_test$y, 
                   x                   = f$den_flat_test$x)
        
  D6 <- data.frame(information_sources = "Prior", 
                   group               = "Test", 
                   y                   = f$den_prior_test$y, 
                   x                   = f$den_prior_test$x)

  D <- as.data.frame(rbind(D4, D5, D6))

  D$information_sources <- factor(D$information_sources, 
                                  levels = (c("Posterior", "Current data", "Prior")))

  post_typeplot <- ggplot(D, aes(x = x, y = y)) +
    geom_line(size = 2, aes(colour = information_sources, lty = information_sources)) +
    theme_bw() + 
    facet_wrap(~group, ncol = 1, scale = "free") + 
    ylab("Density (PDF)") + 
    xlab("values")

      
  densityplot <- ggplot(subset(D, information_sources == "Posterior"), aes(x = x, y = y)) +
    geom_line(size = 2, aes(colour = group)) + 
    ylab("Density (PDF)") + 
    xlab("values") + 
    theme_bw()

  if (two_side == 1) {
      p_value = seq(0, 1, , 100)
      p_value = ifelse(p_value > 0.5, 1 - p_value, p_value)
  }
  if (two_side == 0) {
      p_value = seq(0, 1, , 100)
  }

  Loss_function_test <- pweibull(p_value, 
                                 shape = posterior_test$weibull_shape, 
                                 scale = posterior_test$weibull_scale)*posterior_test$N0_max

  D1 <- data.frame(group = "test", y = Loss_function_test, x = seq(0, 1, , 100))
  D2 <- data.frame(group = c("test"), pvalue = c(posterior_test$pvalue))
  D3 <- data.frame(group = c("test"), pvalue = c(posterior_test$N0_effective))


  lossfun_plot <- ggplot() + 
    geom_line(data = D1, aes(y = y, x = x, colour = group), size = 1) + 
    geom_vline(data = D2, aes(xintercept = pvalue, colour = group), lty = 2) + 
    geom_hline(data = D3, aes(yintercept = pvalue, colour = group), lty = 2) + 
    facet_wrap(~group, ncol = 1) + 
    theme_bw() + 
    ylab("Effective sample size for historical data") + 
    xlab("Bayesian p-value (new vs historical data)")

  if (inequality == "<") {
    hypothesis <- paste("\"We can define mu as the mean for the test", "\n", 
                        "Null Hypothesis (H_0): mu>", H0, "\n", 
                        "Alternative Hypothesis (H_a): mu<", H0, "\n", "\n", 
                        "P(mu<", H0, "|data)=", mean(f$Testpost < H0), "\n", 
                        "We can accept H_a with a Probability of", 
                        mean(f$Testpost < H0))
  }
  
  if (inequality == ">") {
    hypothesis <- paste("\"Define mu as the mean of the data", "\n", 
                        "Null Hypothesis (H_0): mu<", H0, "\n", 
                        "Alternative Hypothesis (H_a): mu>", H0, "\n", "\n", 
                        "P(mu>", H0, "|data)=", mean(f$Testpost > H0), "\n", 
                        "We can accept H_a with a Probability of", 
                        mean(f$Testpost > H0))
  }
    
  ### Print
  prior_for_test_group <- list(`Effective sample size of prior(for test group)` = posterior_test$N0_effective, 
                               `Bayesian p-value (new vs historical data)` = posterior_test$pvalue, 
                               `loss function value` = posterior_test$alpha_loss, 
                               N0_max = posterior_test$N0_max)

  return(list(prior_for_test_group = prior_for_test_group,
              post_typeplot        = post_typeplot,
              densityplot          = densityplot,
              lossfun_plot         = lossfun_plot,
              hypothesis           = hypothesis)
  )

}
