#' @title plot
#' @description normal plot method
#' @import methods
#' @importFrom utils head
#' @importFrom ggplot2 aes_string ggtitle ylim guides guide_legend theme element_blank
#' @importFrom graphics par
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @param x Result
#' @export
setMethod("plot", signature(x = "bdpnormal"), function(x){
  f <- x$f1
  posterior_treatment <- x$posterior_treatment
  posterior_control <- x$posterior_control
  two_side <- x$args1$two_sid
  N0_t <- x$args1$N0_t
  N0_c <- x$args1$N0_c
  arm2 <- x$args1$arm2

  information_sources <- NULL

  if (arm2 == TRUE){
    D1 <- data.frame(information_sources="Posterior",
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

  D4 <- data.frame(information_sources="Posterior",
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

  post_typeplot <- ggplot(D,aes_string(x="x",y="y")) +
    geom_line(size=2,aes_string(colour="information_sources",lty="information_sources")) +
    facet_wrap(~group, ncol=1,scales='free') +
    ylab("Density (PDF)") +
    xlab("Values") +
    theme_bw() +
    ggtitle("Posterior Type Plot")

  densityplot <- ggplot(subset(D,information_sources=="Posterior"),
                        aes_string(x="x",y="y")) +
    geom_line(size=2,aes_string(colour="group")) +
    ylab("Density (PDF)") +
    xlab("Values") +
    theme_bw() +
    ggtitle("Density Plot")


  if(two_side==1){
    p_value <- seq(0,1,length.out=100)
    p_value <- ifelse(p_value>.5,1-p_value,p_value)
  }
  if(two_side==0){
    p_value <- seq(0,1,length.out=100)
  }

  Discount_function_treatment <- pweibull(p_value,
                                          shape=x$args1$weibull_shape[1],
                                          scale=x$args1$weibull_scale[1])
  if(arm2 == TRUE){
    Discount_function_control <- pweibull(p_value,
                                          shape=x$args1$weibull_shape[2],
                                          scale=x$args1$weibull_scale[2])
  }

  D1 <- data.frame(group = "Treatment",
                   y     = Discount_function_treatment,
                   x     = seq(0,1,length.out=100))
  D2 <- data.frame(group=c("Treatment"),pvalue=c(posterior_treatment$pvalue))
  D3 <- data.frame(group=c("Treatment"),pvalue=c(posterior_treatment$alpha_discount))

  if(arm2 == TRUE){
    D4 <- data.frame(group = "Control",
                     y     = Discount_function_control,
                     x     = seq(0,1,length.out=100))
    D5 <- data.frame(group=c("Control"),pvalue=c(posterior_control$pvalue))
    D6 <- data.frame(group=c("Control"),pvalue=c(posterior_control$alpha_discount))
  }


  discountfun_plot <- ggplot()
  if(N0_t!=0){
    discountfun_plot <- discountfun_plot +
      geom_line(data=D1,aes_string(y="y",x="x",colour="group"),size=1) +
      geom_vline(data=D2, aes_string(xintercept="pvalue",colour="group"),lty=2) +
      geom_hline(data=D3, aes_string(yintercept="pvalue",colour="group"),lty=2)
  }
  if(arm2 == TRUE){
    discountfun_plot  <- discountfun_plot +
      geom_line(data=D4,aes_string(y="y",x="x",colour="group"),size=1) +
      geom_vline(data=D5, aes_string(xintercept="pvalue",colour="group"),lty=2) +
      geom_hline(data=D6, aes_string(yintercept="pvalue",colour="group"),lty=2)
  }

  discountfun_plot <- discountfun_plot +
    facet_wrap(~group, ncol=1) +
    theme_bw() +
    ylab("Alpha Discount Value") +
    xlab("Stochastic comparison (New vs Historical Data)") +
    ggtitle("Discount Function") +
    ylim(0,1)

  post_typeplot <- post_typeplot +
    guides(fill=guide_legend(title=NULL)) +
    theme(legend.title=element_blank())

  densityplot <- densityplot +
    guides(fill=guide_legend(title=NULL)) +
    theme(legend.title=element_blank())

  discountfun_plot <- discountfun_plot +
    guides(fill=guide_legend(title=NULL)) +
    theme(legend.title=element_blank())

  op <- par(ask=TRUE)
  plot(post_typeplot)
  plot(densityplot)
  plot(discountfun_plot)
  par(op)
})


#' @title plot
#' @description binomial plot method
#' @import methods
#' @importFrom utils head
#' @importFrom ggplot2 aes_string ggtitle ylim guides guide_legend theme element_blank
#' @importFrom graphics par
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @param x Result
#' @export
setMethod("plot", signature(x = "bdpbinomial"), function(x){
  f                   <- x$f1
  posterior_treatment <- x$posterior_treatment
  posterior_control   <- x$posterior_control
  two_side            <- x$args1$two_side
  N0_t                <- x$args1$N0_t
  N0_c                <- x$args1$N0_c
  N_t                 <- x$args1$N_t
  N_c                 <- x$args1$N_c
  arm2                <- x$args1$arm2

  information_sources <- NULL

  ### Create dataframes for plotting via ggplot2
  D1 <- D2 <- D3 <- D5 <- D6 <- NULL

  if(arm2){
    D1 <- data.frame(information_sources='Posterior',
                     group="Control",
                     y=f$density_post_control$y,
                     x=f$density_post_control$x)

    if(!is.null(f$density_flat_control)){
      D2 <- data.frame(information_sources="Current Data",
                       group="Control",
                       y=f$density_flat_control$y,
                       x=f$density_flat_control$x)
    }

    if(!is.null(f$density_prior_control)){
      D3 <- data.frame(information_sources="Prior",
                       group="Control",
                       y=f$density_prior_control$y,
                       x=f$density_prior_control$x)
    }
  }

  D4 <- data.frame(information_sources='Posterior',
                   group="Treatment",
                   y=f$density_post_treatment$y,
                   x=f$density_post_treatment$x)

  if(!is.null(f$density_flat_treatment)){
    D5 <- data.frame(information_sources="Current Data",
                     group="Treatment",
                     y=f$density_flat_treatment$y,
                     x=f$density_flat_treatment$x)
  }

  if(!is.null(f$density_prior_treatment)){
    D6 <- data.frame(information_sources="Prior",
                     group="Treatment",
                     y=f$density_prior_treatment$y,
                     x=f$density_prior_treatment$x)
  }

  D <- as.data.frame(rbind(D1,D2,D3,D4,D5,D6))

  D$information_sources <- factor(D$information_sources,
                                  levels = (c("Posterior","Current Data","Prior")))

  ##############################################################################
  ### Posterior Type Plots
  ##############################################################################
  post_typeplot <- ggplot(D,aes_string(x="x",y="y")) +
    geom_line(size=2,aes_string(color="information_sources",lty="information_sources")) +
    theme_bw() +
    facet_wrap(~group, ncol=1,scales='free') +
    ylab("Density (PDF)") +
    xlab("Values") +
    ggtitle("Posterior Type Plot") +
    guides(fill=guide_legend(title=NULL)) +
    theme(legend.title=element_blank())

  ##############################################################################
  ### Density Plots
  ##############################################################################
  densityplot <- ggplot(subset(D,information_sources=="Posterior"),
                        aes_string(x="x",y="y")) +
    geom_line(size=2,aes_string(color="group")) +
    ylab("Density (PDF)") +
    xlab("Values") +
    theme_bw() +
    ggtitle("Density Plot") +
    guides(fill=guide_legend(title=NULL)) +
    theme(legend.title=element_blank())


  ##############################################################################
  ### Discount function plot
  ### - Only makes sense to plot if Current/historical treatment are present or
  ###   both current/historical control are present
  ##############################################################################
  if(two_side){
    p_value <- seq(0,1,length.out=100)
    p_value <- ifelse(p_value>.5,1-p_value,p_value)
  } else{
    p_value <- seq(0,1,length.out=100)
  }


  discountfun_plot <- NULL

  if(!is.null(N0_t) & !is.null(N_t)){
    discount_function_treatment <- pweibull(p_value,
                                            shape=x$args1$weibull_shape[1],
                                            scale=x$args1$weibull_scale[1])
    D1 <- data.frame(group = "Treatment",
                     y     = discount_function_treatment,
                     x     = seq(0,1,length.out=100))
    D2 <- data.frame(group="Treatment", pvalue=c(posterior_treatment$pvalue))
    D3 <- data.frame(group="Treatment", pvalue=c(posterior_treatment$alpha_discount))

    discountfun_plot <- ggplot() +
      geom_line(data=D1,aes_string(y="y",x="x",color="group"),size=1) +
      geom_vline(data=D2, aes_string(xintercept="pvalue",color="group"),lty=2) +
      geom_hline(data=D3, aes_string(yintercept="pvalue",color="group"),lty=2)
  }


  if(arm2){
    if(!is.null(N0_c) & !is.null(N_c)){
      discount_function_control <- pweibull(p_value,
                                            shape=x$args1$weibull_shape[2],
                                            scale=x$args1$weibull_scale[2])

      D4 <- data.frame(group = "Control",
                       y     = discount_function_control,
                       x     = seq(0,1,length.out=100))
      D5 <- data.frame(group="Control",pvalue=c(posterior_control$pvalue))
      D6 <- data.frame(group="Control",pvalue=c(posterior_control$alpha_discount))

      discountfun_plot  <- discountfun_plot +
        geom_line(data=D4,aes_string(y="y",x="x",color="group"),size=1) +
        geom_vline(data=D5,aes_string(xintercept="pvalue",color="group"),lty=2) +
        geom_hline(data=D6,aes_string(yintercept="pvalue",color="group"),lty=2)
    }
  }


  if(!is.null(discountfun_plot)){
    discountfun_plot <- discountfun_plot +
      facet_wrap(~group, ncol=1) +
      theme_bw() +
      ylab("Alpha Discount Value") +
      xlab("Stochastic comparison (New vs Historical Data)") +
      ggtitle("Discount Function") +
      ylim(0,1) +
      guides(fill=guide_legend(title=NULL)) +
      theme(legend.title=element_blank())
  }


  op <- par(ask=TRUE)
  plot(post_typeplot)
  plot(densityplot)
  if(!is.null(discountfun_plot)){
    plot(discountfun_plot)
  }
  par(op)
})


#' @title plot
#' @description survival plot method
#' @import methods
#' @importFrom utils head
#' @importFrom ggplot2 aes_string ggtitle ylim guides guide_legend theme element_blank
#' @importFrom graphics par
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @param x Result
#' @export
setMethod("plot", signature(x = "bdpsurvival"), function(x){
  f                   <- x$f1
  posterior_treatment <- x$posterior_treatment
  two_side            <- x$args1$two_side

  information_sources <- NULL


  D4 <- data.frame(information_sources = "Posterior",
                   group               = "Treatment",
                   y                   = f$density_post_treatment$y,
                   x                   = f$density_post_treatment$x)

  D5 <- data.frame(information_sources = "Current Data",
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

  post_typeplot <- ggplot(D, aes_string(x = "x", y = "y")) +
    geom_line(size = 2, aes_string(colour = "information_sources", lty = "information_sources")) +
    theme_bw() +
    facet_wrap(~group, ncol = 1, scales = "free") +
    ylab("Density (PDF)") +
    xlab("Values") +
    ggtitle("Posterior Type Plot")

  densityplot <- ggplot(subset(D, information_sources == "Posterior"), aes_string(x = "x", y = "y")) +
    geom_line(size = 2, aes_string(colour = "group")) +
    ylab("Density (PDF)") +
    xlab("Values") +
    theme_bw() +
    ggtitle("Density Plot")


  if (two_side == 1) {
    p_value = seq(0, 1, length.out=100)
    p_value = ifelse(p_value > 0.5, 1 - p_value, p_value)
  }
  if (two_side == 0) {
    p_value = seq(0, 1, length.out=100)
  }

  Loss_function_treatment <- pweibull(p_value,
                                      shape = x$args1$weibull_shape[1],
                                      scale = x$args1$weibull_scale[1])

  D1 <- data.frame(group = "Treatment",
                   y     = Loss_function_treatment,
                   x     = seq(0, 1, length.out=100))
  D2 <- data.frame(group = "Treatment", pvalue = posterior_treatment$pvalue)
  D3 <- data.frame(group = "Treatment", pvalue = posterior_treatment$alpha_discount)


  discountfun_plot <- ggplot() +
    geom_line(data = D1, aes_string(y = "y", x = "x", colour = "group"), size = 1) +
    geom_vline(data = D2, aes_string(xintercept = "pvalue", colour = "group"), lty = 2) +
    geom_hline(data=D3, aes_string(yintercept ="pvalue", colour="group"),lty=2)

  discountfun_plot <- discountfun_plot +
    facet_wrap(~group, ncol = 1) +
    theme_bw() +
    ylab("Alpha Discount Value") +
    xlab("Stochastic comparison (New vs Historical Data)") +
    ggtitle("Discount Function") +
    ylim(0,1)

  post_typeplot <- post_typeplot +
    guides(fill=guide_legend(title=NULL)) +
    theme(legend.title=element_blank())

  densityplot <- densityplot +
    guides(fill=guide_legend(title=NULL)) +
    theme(legend.title=element_blank())

  discountfun_plot <- discountfun_plot +
    guides(fill=guide_legend(title=NULL)) +
    theme(legend.title=element_blank())

  op <- par(ask=TRUE)
  plot(post_typeplot)
  plot(densityplot)
  plot(discountfun_plot)
  par(op)
})


