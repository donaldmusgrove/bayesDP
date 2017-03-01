#' plot
#' plot
#' @title plot: plot
#' @param x Results
#' @importFrom graphics par
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @rdname plot
#' @aliases plot, bdpnormal, bdpbinomial, bdpregression_linear, bdpsurvival
#' @export plot
setMethod("plot", signature(x = "bdpnormal"), function(x){
  f <- x$f1
  posterior_treatment <- x$posterior_treatment
  posterior_control <- x$posterior_control
  two_side <- x$args1$two_sid
  N0_t <- x$args1$N0_t
  N0_c <- x$args1$N0_c
  arm2 <- x$args1$arm2
  if (arm2 == TRUE){
    D1 <- data.frame(information_sources='Posterior',
                     group="Control",
                     y=f$den_post_control$y,
                     x=f$den_post_control$x)
    D2 <- data.frame(information_sources="Current Data",
                     group="Control",
                     y=f$den_flat_control$y,
                     x=f$den_flat_control$x)
    D3 <- data.frame(information_sources="Prior",
                     group="Control",
                     y=f$den_prior_control$y,
                     x=f$den_prior_control$x)
  }

  D4 <- data.frame(information_sources='Posterior',
                   group="Test",
                   y=f$den_post_treatment$y,
                   x=f$den_post_treatment$x)
  D5 <- data.frame(information_sources="Current Data",
                   group="Test",
                   y=f$den_flat_treatment$y,
                   x=f$den_flat_treatment$x)
  D6 <- data.frame(information_sources="Prior",
                   group="Test",
                   y=f$den_prior_treatment$y,
                   x=f$den_prior_treatment$x)

  if(is.null(N0_t) == TRUE & is.null(N0_c) == TRUE){
    D <- as.data.frame(rbind(D4,D5))
  }
  if(is.null(N0_t) == TRUE & is.null(N0_c) == FALSE){
    D <- as.data.frame(rbind(D4,D5,D1,D2,D3))
  }
  if(is.null(N0_t) == FALSE & is.null(N0_c) == TRUE){
    D <- as.data.frame(rbind(D4,D5,D6))
  }
  if(is.null(N0_t) == FALSE & is.null(N0_c) == FALSE){
    D <- as.data.frame(rbind(D4,D5,D6,D1,D2,D3))
  }

  D$information_sources <- factor(D$information_sources,
                                  levels = (c("Posterior","Current Data","Prior")))

  post_typeplot <- ggplot(D,aes(x=x,y=y)) +
    geom_line(size=2,aes(colour=information_sources,lty=information_sources)) +
    facet_wrap(~group, ncol=1,scales='free') +
    ylab("Density (PDF)") +
    xlab("Values") +
    theme_bw() +
    ggtitle("Posterior Type Plot")

  densityplot <- ggplot(subset(D,information_sources=="Posterior"),
                        aes(x=x,y=y)) +
    geom_line(size=2,aes(colour=group)) +
    ylab("Density (PDF)") +
    xlab("Values") +
    theme_bw() +
    ggtitle("Density Plot")


  if(two_side==1){
    p_value <- seq(0,1,,100)
    p_value <- ifelse(p_value>.5,1-p_value,p_value)
  }
  if(two_side==0){
    p_value <- seq(0,1,,100)
  }

  Discount_function_treatment <- pweibull(p_value,
                                          shape=posterior_treatment$weibull_shape,
                                          scale=posterior_treatment$weibull_scale)
  if(arm2 == TRUE){
    Discount_function_control <- pweibull(p_value,
                                          shape=posterior_control$weibull_shape,
                                          scale=posterior_control$weibull_scale)
  }

  D1 <- data.frame(group="Treatment",y=Discount_function_treatment,x=seq(0,1,,100))
  D2 <- data.frame(group=c("Treatment"),pvalue=c(posterior_treatment$pvalue))
  D3 <- data.frame(group=c("Treatment"),pvalue=c(posterior_treatment$alpha_discount))

  if(arm2 == TRUE){
    D4 <- data.frame(group="Control",y=Discount_function_control,x=seq(0,1,,100))
    D5 <- data.frame(group=c("Control"),pvalue=c(posterior_control$pvalue))
    D6 <- data.frame(group=c("Control"),pvalue=c(posterior_control$alpha_discount))
  }


  discountfun_plot <- ggplot()
  if(N0_t!=0){
    discountfun_plot <- discountfun_plot +
      geom_line(data=D1,aes(y=y,x=x,colour=group),size=1) +
      geom_vline(data=D2, aes(xintercept =pvalue,colour=group),lty=2) +
      geom_hline(data=D3, aes(yintercept =pvalue,colour=group),lty=2)
  }
  if(arm2 == TRUE){
    discountfun_plot  <- discountfun_plot +
      geom_line(data=D4,aes(y=y,x=x,colour=group),size=1) +
      geom_vline(data=D5, aes(xintercept =pvalue,colour=group),lty=2) +
      geom_hline(data=D6, aes(yintercept =pvalue,colour=group),lty=2)
  }

  discountfun_plot <- discountfun_plot +
    facet_wrap(~group, ncol=1) +
    theme_bw() +
    ylab("Alpha Discount Value") +
    xlab("Bayesian p-value (New vs Historical Data)") +
    ggtitle("Discount Function") +
    ylim(0,1)

  post_typeplot <- post_typeplot + guides(fill=guide_legend(title=NULL))

  post_typeplot <- post_typeplot + theme(legend.title=element_blank())

  densityplot <- densityplot + guides(fill=guide_legend(title=NULL))

  densityplot <- densityplot + theme(legend.title=element_blank())

  discountfun_plot <- discountfun_plot + guides(fill=guide_legend(title=NULL))

  discountfun_plot <- discountfun_plot + theme(legend.title=element_blank())

  op <- par(ask=TRUE)
  plot(post_typeplot)
  plot(densityplot)
  plot(discountfun_plot)
  par(op)
})


setMethod("plot", signature(x = "bdpbinomial"), function(x){
  f                   <- x$f1
  posterior_treatment <- x$posterior_treatment
  posterior_control   <- x$posterior_control
  two_side            <- x$args1$two_side
  N0_t                <- x$args1$N0_t
  N0_c                <- x$args1$N0_c
  arm2                <- x$args1$arm2
  if (arm2){
    D1 <- data.frame(information_sources='Posterior',
                     group="Control",
                     y=f$density_post_control$y,
                     x=f$density_post_control$x)

    D2 <- data.frame(information_sources="Current Data",
                     group="Control",
                     y=f$density_flat_control$y,
                     x=f$density_flat_control$x)

    D3 <- data.frame(information_sources="Prior",
                     group="Control",
                     y=f$density_prior_control$y,
                     x=f$density_prior_control$x)
  }

  D4 <- data.frame(information_sources='Posterior',
                   group="Treatment",
                   y=f$density_post_treatment$y,
                   x=f$density_post_treatment$x)

  D5 <- data.frame(information_sources="Current Data",
                   group="Treatment",
                   y=f$density_flat_treatment$y,
                   x=f$density_flat_treatment$x)

  D6 <- data.frame(information_sources="Prior",
                   group="Treatment",
                   y=f$density_prior_treatment$y,
                   x=f$density_prior_treatment$x)

  if(is.null(N0_t) == TRUE & is.null(N0_c) == TRUE){
    D <- as.data.frame(rbind(D4,D5))
  }
  if(is.null(N0_t) == TRUE & is.null(N0_c) == FALSE){
    D <- as.data.frame(rbind(D4,D5,D1,D2,D3))
  }
  if(is.null(N0_t) == FALSE & is.null(N0_c) == TRUE){
    D <- as.data.frame(rbind(D4,D5,D6))
  }
  if(is.null(N0_t) == FALSE & is.null(N0_c) == FALSE){
    D <- as.data.frame(rbind(D4,D5,D6,D1,D2,D3))
  }

  D$information_sources <- factor(D$information_sources,
                                  levels = (c("Posterior","Current Data","Prior")))

  post_typeplot <- ggplot(D,aes(x=x,y=y)) +
    geom_line(size=2,aes(color=information_sources,lty=information_sources)) +
    theme_bw() +
    facet_wrap(~group, ncol=1,scales='free') +
    ylab("Density (PDF)") +
    xlab("Values") +
    ggtitle("Posterior Type Plot")

  densityplot <- ggplot(subset(D,information_sources=="Posterior"),
                        aes(x=x,y=y)) +
    geom_line(size=2,aes(color=group)) +
    ylab("Density (PDF)") +
    xlab("Values") +
    theme_bw() +
    ggtitle("Density Plot")


  if(two_side==1){
    p_value <- seq(0,1,,100)
    p_value <- ifelse(p_value>.5,1-p_value,p_value)
  }
  if(two_side==0){
    p_value <- seq(0,1,,100)
  }

  discount_function_treatment <- pweibull(p_value,
                                          shape=posterior_treatment$weibull_shape,
                                          scale=posterior_treatment$weibull_scale)
  if(arm2 == TRUE){
    discount_function_control <- pweibull(p_value,
                                          shape=posterior_control$weibull_shape,
                                          scale=posterior_control$weibull_scale)
  }

  D1 <- data.frame(group="Treatment",y=discount_function_treatment,x=seq(0,1,,100))
  D2 <- data.frame(group="Treatment",pvalue=c(posterior_treatment$pvalue))
  D3 <- data.frame(group="Treatment",pvalue=c(posterior_treatment$alpha_discount))
  if (arm2 == TRUE){
    D4 <- data.frame(group="Control",y=discount_function_control,x=seq(0,1,,100))
    D5 <- data.frame(group="Control",pvalue=c(posterior_control$pvalue))
    D6 <- data.frame(group="Control",pvalue=c(posterior_control$alpha_discount))
  }


  discountfun_plot <- ggplot()
  if(N0_t!=0){
    discountfun_plot <- discountfun_plot +
      geom_line(data=D1,aes(y=y,x=x,color=group),size=1) +
      geom_vline(data=D2, aes(xintercept=pvalue,color=group),lty=2) +
      geom_hline(data=D3, aes(yintercept=pvalue,color=group),lty=2)
  }
  if(arm2 == TRUE){
    discountfun_plot  <- discountfun_plot +
      geom_line(data=D4,aes(y=y,x=x,color=group),size=1) +
      geom_vline(data=D5,aes(xintercept=pvalue,color=group),lty=2) +
      geom_hline(data=D6,aes(yintercept=pvalue,color=group),lty=2)
  }

  discountfun_plot <- discountfun_plot +
    facet_wrap(~group, ncol=1) +
    theme_bw() +
    ylab("Alpha Discount Value") +
    xlab("Bayesian p-value (New vs Historical Data)") +
    ggtitle("Discount Function") +
    ylim(0,1)

  post_typeplot <- post_typeplot + guides(fill=guide_legend(title=NULL))

  post_typeplot <- post_typeplot + theme(legend.title=element_blank())

  densityplot <- densityplot + guides(fill=guide_legend(title=NULL))

  densityplot <- densityplot + theme(legend.title=element_blank())

  discountfun_plot <- discountfun_plot + guides(fill=guide_legend(title=NULL))

  discountfun_plot <- discountfun_plot + theme(legend.title=element_blank())

  op <- par(ask=TRUE)
  plot(post_typeplot)
  plot(densityplot)
  plot(discountfun_plot)
  par(op)
})


setMethod("plot", signature(x = "bdpregression_linear"), function(x){
  f          <- x$f1
  posterior  <- x$est
  two_side   <- x$args1$two_side

  D4 <- data.frame(information_sources = "Posterior",
                   y                   = f$den_post$y,
                   x                   = f$den_post$x)

  D5 <- data.frame(information_sources = "Current Data",
                   y                   = f$den_flat$y,
                   x                   = f$den_flat$x)

  D6 <- data.frame(information_sources = "Prior",
                   y                   = f$den_prior$y,
                   x                   = f$den_prior$x)

  D <- as.data.frame(rbind(D4, D5, D6))

  D$information_sources <- factor(D$information_sources,
                                  levels = (c("Posterior", "Current Data", "Prior")))

  post_typeplot <- ggplot(D, aes(x = x, y = y)) +
    geom_line(size = 2, aes(colour = information_sources, lty = information_sources)) +
    theme_bw() +
    ylab("Density (PDF)") +
    xlab("Values") +
    ggtitle("Posterior Type Plot")


  densityplot <- ggplot(subset(D, information_sources == "Posterior"), aes(x = x, y = y)) +
    geom_line(size = 2) +
    ylab("Density (PDF)") +
    xlab("Values") +
    theme_bw() +
    ggtitle("Density Plot")

  if (two_side == 1) {
    p_value = seq(0, 1, , 100)
    p_value = ifelse(p_value > 0.5, 1 - p_value, p_value)
  }
  if (two_side == 0) {
    p_value = seq(0, 1, , 100)
  }

  Loss_function <- pweibull(p_value,
                            shape = posterior$weibull_shape,
                            scale = posterior$weibull_scale)

  D1 <- data.frame(y = Loss_function, x = seq(0, 1, , 100))
  D2 <- data.frame(pvalue = c(posterior$pvalue))

  lossfun_plot <- ggplot() +
    geom_line(data = D1, aes(y = y, x = x), size = 1) +
    geom_vline(data = D2, aes(xintercept = pvalue), lty = 2) +
    theme_bw() +
    ylab("Effective Sample Size for Historical Data") +
    xlab("Bayesian p-value (New vs Historical Data)") +
    ggtitle("Discount Function Plot")

  post_typeplot <- post_typeplot + guides(fill=guide_legend(title=NULL))

  post_typeplot <- post_typeplot + theme(legend.title=element_blank())

  densityplot <- densityplot + guides(fill=guide_legend(title=NULL))

  densityplot <- densityplot + theme(legend.title=element_blank())

  discountfun_plot <- discountfun_plot + guides(fill=guide_legend(title=NULL))

  discountfun_plot <- discountfun_plot + theme(legend.title=element_blank())

  op <- par(ask=TRUE)
  plot(post_typeplot)
  plot(densityplot)
  plot(lossfun_plot)
  par(op)
})


setMethod("plot", signature(x = "bdpsurvival"), function(x){
  f                   <- x$f1
  posterior_treatment <- x$posterior_treatment
  two_side            <- x$args1$two_side
  starts              <- c(0,x$args1$breaks)

  D4 <- data.frame(information_sources = "Posterior",
                   group               = "Treatment",
                   y                   = f$density_post_treatment$y,
                   x                   = f$density_post_treatment$x)

  D5 <- data.frame(information_sources = "Current data",
                   group               = "Treatment",
                   y                   = f$density_flat_treatment$y,
                   x                   = f$density_flat_treatment$x)

  D6 <- data.frame(information_sources = "Prior",
                   group               = "Treatment",
                   y                   = f$density_prior_treatment$y,
                   x                   = f$density_prior_treatment$x)

  D <- as.data.frame(rbind(D4, D5, D6))

  D$information_sources <- factor(D$information_sources,
                                  levels = (c("Posterior", "Current Data", "Prior")))

  D$start <- paste0("Interval start: ", D$start)

  post_typeplot <- ggplot(D, aes(x = x, y = y)) +
    geom_line(size = 2, aes(colour = information_sources, lty = information_sources)) +
    theme_bw() +
    facet_wrap(~group+start, ncol = 1, scales = "free") +
    ylab("Density (PDF)") +
    xlab("Values") +
    ggtitle("Posterior Type Plot")

  densityplot <- ggplot(subset(D, information_sources == "Posterior"), aes(x = x, y = y)) +
    geom_line(size = 2, aes(colour = group)) +
    ylab("Density (PDF)") +
    xlab("Values") +
    theme_bw() +
    ggtitle("Density Plot")


  if (two_side == 1) {
    p_value = seq(0, 1, , 100)
    p_value = ifelse(p_value > 0.5, 1 - p_value, p_value)
  }
  if (two_side == 0) {
    p_value = seq(0, 1, , 100)
  }

  Loss_function_treatment <- pweibull(p_value,
                                      shape = posterior_treatment$weibull_shape,
                                      scale = posterior_treatment$weibull_scale)

  D1 <- data.frame(group = "treatment", y = Loss_function_treatment, x = seq(0, 1, , 100))
  D2 <- data.frame(group = c("treatment"), pvalue = c(posterior_treatment$pvalue))

  lossfun_plot <- ggplot() +
    geom_line(data = D1, aes(y = y, x = x, colour = group), size = 1) +
    geom_vline(data = D2, aes(xintercept = pvalue, colour = group), lty = 2) +
    facet_wrap(~group, ncol = 1) +
    theme_bw() +
    ylab("Effective Sample Size for Historical Data") +
    xlab("Bayesian p-value (New vs Historical Data)") +
    ggtitle("Discount Function Plot")

  post_typeplot <- post_typeplot + guides(fill=guide_legend(title=NULL))

  post_typeplot <- post_typeplot + theme(legend.title=element_blank())

  densityplot <- densityplot + guides(fill=guide_legend(title=NULL))

  densityplot <- densityplot + theme(legend.title=element_blank())

  lossfun_plot <- lossfun_plot + guides(fill=guide_legend(title=NULL))

  lossfun_plot <- lossfun_plot + theme(legend.title=element_blank())

  op <- par(ask=TRUE)
  plot(post_typeplot)
  plot(densityplot)
  plot(lossfun_plot)
  par(op)
})


