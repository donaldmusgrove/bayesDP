#normalrtc(test,control)


#loss_function_inputs,test_group_info,control_group_info,number_mcmc

loss_function_inputs = list(test_group_info = list(N0_max = 10,
                                                 weibull_scale = 0.1,
                                                 weibull_shape = 2),
                            control_group_info = list(N0_max = 10,
                                                      weibull_scale = 0.1,
                                                      weibull_shape = 2))
test_group_info = list(test_group_info = list(mu = 10,
                                              sigma = 0.1,
                                              N = 2))
                       control_group_info = list(mu = 10,
                                                 sigma = 0.1,
                                                 N = 2),
control_group_info = list(test_group_info = list(mu = 10,
                                                 sigma = 0.1,
                                                 N = 2),
                          control_group_info = list(mu = 10,
                                                    sigma = 0.1,
                                                    N = 2))
two_sided_loss_function = False,
number_mcmc=1000
#c(10,0.1,2),c(10,0.1,2)c(10,1,10),c(10,1,10)c(10,1,10),c(10,1,10)
