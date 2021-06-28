gg_target <- function(mods, obs, colorval = NULL, colorval.name = NULL, axis_begin, axis_end, by, label){
  
  # compute normalized ME, Eq. 9 Jolliff et al., (2009) - Y-axis
  model_nME <- plyr::laply(mods, .fun = function(x, y){(mean(y) - mean(x))/sd(y)}, y = obs)
  
  # compute the standard deviation difference (sigma_d), Eq. 11 Jolliff et al., (2009) - X-axis first part
  model_sigmaD <- plyr::laply(mods, .fun = function(x, y){sign(sd(x) - sd(y))}, y = obs)
  
  # if sign not negative then assign 1 to zero
  model_sigmaD[model_sigmaD == 0] <- 1
  
  # compute the normalized unbiased RMSD, called (model_uRMSDnorm below) RMSD*' in Eq. 8 Jolliff et al., (2009) - X-axis second part
  model_uRMSDnorm <- plyr::laply(mods, .fun = function(x, y){   
    # Eq. 2 Jolliff
    sigma_ast <- sd(x)/sd(y)
    corxy <- cor(x,y)
    # account for possible sd = 0 (mean field), then cor = 0
    corxy[is.na(corxy)] <- 0
    uRMSDnorm <- sqrt(1 + sigma_ast^2 - 2*sigma_ast*corxy)}, 
    y = obs)
  
  # additional values to be plotted as colour on the diagram, the r
  model_cors <- plyr::laply(mods, .fun = function(x,y){cor(x,y)}, y = obs)
  # account for possible sd=0 (mean field), then cor = 0
  model_cors[is.na(model_cors)] <- 0
  
  if(is.null(colorval)){colorval = model_cors}
  if(is.null(colorval.name)){colorval.name = 'Correlation'}
  
  
  data <- data.frame(
    nME = model_nME,
    sigmaD = model_sigmaD,
    uRMSDnorm = model_uRMSDnorm,
    uRMSDnorm_sigmaD = model_uRMSDnorm*model_sigmaD,
    r = model_cors, 
    colvar = colorval,
    model = factor(names(mods))
  )
  
  # make circle at some cuts
  # the cuts are determined based on the correlation, see Jollief et al., (2009), Eq. 14
  cuts = c(0.44, 0.71) 
  circle <- expand.grid(theta = seq(0, 2 * pi, length = 100), 
                        r = cuts)
  circle$x <- with(circle, r * sin(theta))
  circle$y <- with(circle, r * cos(theta))
  
  # make circle at some cuts
  cuts = c(1)
  circle2 <- expand.grid(theta = seq(0, 2 * pi, length = 100), 
                         r = cuts)
  circle2$x <- with(circle2, r * sin(theta))
  circle2$y <- with(circle2, r * cos(theta))
  
  # labels of the semicircles
  radius <- unique(c(circle$r, unique(circle2$r)))
  labels <- data.frame(x = 0, y = -radius,
                       lbs = c(0.9, 0.7, 1))
  labels$x[1] <- circle[circle$r == 0.44,]$x[37]
  labels$y[1] <- circle[circle$r == 0.44,]$y[37]
  labels$x[2] <- circle[circle$r == 0.71,]$x[37]
  labels$y[2] <- circle[circle$r == 0.71,]$y[37]
  labels$x[3] <- circle2[circle2$r == 1,]$x[37]
  labels$y[3] <- circle2[circle2$r == 1,]$y[37]
  
  ### prepare the x and y axis
  # constants
  axis_begin  = axis_begin
  axis_end    = axis_end
  by = by
  
  # point to plot
  my_point <- data.frame(x=1, y=1)
  
  # chart 
  tick_frame <- data.frame(ticks = seq(axis_begin, axis_end, by = by), 
                           zero=0) %>%
    subset(ticks != 0)
  
  lab_frame <- data.frame(lab = range(axis_begin, axis_end),zero = 0) %>%
    subset(lab != 0)
  
  tick_sz <- (tail(lab_frame$lab, 1) -  lab_frame$lab[1]) / 128
  
  # make the ggplot
  ggplot(data = data,  aes(x = uRMSDnorm_sigmaD, y = nME)) +
    geom_path(aes_string(x = 'x', y = 'y', fill = NULL),
              data = circle, col = 'black',linetype = 2) +
    geom_path(aes_string(x = 'x', y = 'y', fill = NULL),
              data = circle2, col = 'black', size = 0.8)+
    geom_text(aes_string(x = 'x', y = 'y', label = 'lbs',
                         vjust = 1, hjust = -0.1, fill = NULL),
              size = 4, data = labels, col = 'black') +
    geom_segment(x = 0, xend = 0, 
                 y = lab_frame$lab[1], yend = tail(lab_frame$lab, 1),
                 size = 0.5) +
    # x axis line
    geom_segment(y = 0, yend = 0, 
                 x = lab_frame$lab[1], xend = tail(lab_frame$lab, 1),
                 size = 0.5) +
    # x ticks
    geom_segment(data = tick_frame, 
                 aes(x = ticks, xend = ticks, 
                     y = zero, yend = zero + tick_sz)) +
    # y ticks
    geom_segment(data = tick_frame, 
                 aes(x = zero, xend = zero + tick_sz, 
                     y = ticks, yend = ticks)) + 
    # labels
    geom_text(data=lab_frame, aes(x=lab, y=zero, label=lab),vjust=1.5, size = 4) +
    geom_text(data=lab_frame, aes(x=zero, y=lab, label=lab),hjust=1.5, size = 4) + 
    geom_point(aes(fill = colvar), pch = 21, size = 7) +
    ylab( TeX('$SDE^* \\cdot \\sign(\\sigma_d)') ) + 
    xlab(TeX('ME^*')) +                  
    {if(label == TRUE)geom_label_repel(aes(label = model),
                                       box.padding   = 0.35, 
                                       point.padding = 0.5,
                                       segment.color = 'grey50', 
                                       size = 5)} + 
    theme_classic() +
    theme(plot.margin=unit(c(1,1,1,1),"cm"))+
    labs(fill=colorval.name) +
    theme(axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          text=element_text(size=16,  family="Palatino"),
          axis.title.x = element_text(margin = margin(t = -17, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = -17, b = 0, l = 0)),
          legend.position="right",
          legend.justification="right",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-20), 
          legend.text.align = 1, 
          legend.title = element_text(vjust = 3))  + 
    scale_fill_viridis(option = 'A')
  
}

