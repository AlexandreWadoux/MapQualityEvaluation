gg_solar <- function(mods, obs, colorval = NULL, colorval.name = NULL, x.axis_begin, x.axis_end, y.axis_end, by, label){
  
  # compute standardized ME
  model_nME <- plyr::laply(mods, .fun = function(x, y){(mean(y) - mean(x))/sd(y)}, y = obs)
  
  # compute the standardized SDE, called (model_uRMSDnorm below) SDE*
  model_uRMSDnorm <- plyr::laply(mods, .fun = function(x, y){   
    # Eq. 2 Jolliff
    sigma_ast <- sd(x)/sd(y)
    corxy <- cor(x,y)
    # account for possible sd = 0 (mean field), then cor = 0
    corxy[is.na(corxy)] <- 0
    uRMSDnorm <- sqrt(1 + sigma_ast^2 - 2*sigma_ast*corxy)}, 
    y = obs)
  
  # correlation r
  model_cors <- plyr::laply(mods, .fun = function(x,y){cor(x,y)}, y = obs)
  # account for possible sd = 0 (mean field), then cor = 0
  model_cors[is.na(model_cors)] <- 0
  
  if(is.null(colorval)){colorval = model_cors}
  if(is.null(colorval.name)){colorval.name = 'Correlation'}
  
  data <- data.frame(
    nME = model_nME,
    sigmaD = 1,
    uRMSDnorm = model_uRMSDnorm,
    uRMSDnorm_sigmaD = model_uRMSDnorm,
    r = model_cors, 
    colvar = colorval,
    model = factor(names(mods))
  )
  
  # make circle at some cuts
  # the cuts are determined based on the correlation, see Jollief et al., (2009), Eq. 14
  cuts = c(0.31, 0.44, 0.71)
  circle <- plyr::adply(cuts, 1, .fun = function(x){data.frame(x = x * cos(seq(0,pi,pi/1000)),
                                                               y = x * sin(seq(0,pi,pi/1000)),
                                                               label = x)})
  
  
  # make circle at some cuts
  cuts = c(1)
  circle2 <- plyr::adply(cuts, 1, .fun = function(x){data.frame(x = x * cos(seq(0,pi,pi/1000)),
                                                                y = x * sin(seq(0,pi,pi/1000)),
                                                                label = x)})
  
  ### prepare the x and y axis
  # point to plot
  my_point <- data.frame(x=1, y=1)
  
  # chart junk data
  x.tick_frame <- data.frame(ticks = seq(x.axis_begin, x.axis_end, by = by), 
                             zero=0)
  y.tick_frame <- data.frame(ticks = seq(0, y.axis_end, by = by), 
                             zero=0)
  x.lab_frame <- data.frame(lab = c(0, range(x.axis_begin, x.axis_end)),zero = 0)
  y.lab_frame <- data.frame(lab = range(0, y.axis_end),zero = 0)
  
  tick_sz <- (tail(x.lab_frame$lab, 1) -  x.lab_frame$lab[1]) / 128
  
  # make polygons for yellow color 
  cuts = c(0, 0.31)
  circle095 <- plyr::adply(cuts, 1, .fun = function(x){data.frame(x = x * cos(seq(0,pi,pi/1000)),
                                                                  y = x * sin(seq(0,pi,pi/1000)),
                                                                  label = x)})
  cuts = c(0.31, 0.44)
  circle09 <- plyr::adply(cuts, 1, .fun = function(x){data.frame(x = x * cos(seq(0,pi,pi/1000)),
                                                                 y = x * sin(seq(0,pi,pi/1000)),
                                                                 label = x)})
  cuts = c(0.44, 0.71)
  circle07 <- plyr::adply(cuts, 1, .fun = function(x){data.frame(x = x * cos(seq(0,pi,pi/1000)),
                                                                 y = x * sin(seq(0,pi,pi/1000)),
                                                                 label = x)})
  cuts = c(0.71, 1)
  circle0 <- plyr::adply(cuts, 1, .fun = function(x){data.frame(x = x * cos(seq(0,pi,pi/1000)),
                                                                y = x * sin(seq(0,pi,pi/1000)),
                                                                label = x)})
  circle0$lab <- as.factor('r>0')
  circle07$lab <- as.factor('r>0.7')
  circle09$lab <- as.factor('r>0.9')
  circle095$lab <- as.factor('r>0.95')
  
  circ_pol <- rbind(circle0, circle07, circle09, circle095)
  
  # make the ggplot
  base.plot <- ggplot(data = data,  aes(x = nME, y = uRMSDnorm_sigmaD)) +
    geom_polygon(data = circ_pol, aes(x, y, fill = lab)) + scale_fill_manual(values = c('r>0' = '#FEECA4','r>0.7' = '#FEF4B6', 'r>0.9' = '#FFFCD7', 'r>0.95'=  '#FFFFE5'),
                                                                             name = '',
                                                                             labels = expression(italic(r)>0,italic(r)>0.7, italic(r)>0.9, italic(r)>0.95))+
    geom_path(data = circle2, aes(x = x, y = y, group = label), colour = "black", size = 0.8) +
    geom_segment(x = 0, xend = 0, 
                 y = 0, yend = tail(y.lab_frame$lab, 1),
                 size = 0.5) +
    # x axis line
    geom_segment(y = 0, yend = 0, 
                 x = x.axis_begin, xend = x.axis_end,
                 size = 0.5) +
    # x ticks
    geom_segment(data = x.tick_frame, 
                 aes(x = ticks, xend = ticks, 
                     y = zero, yend = zero + tick_sz)) +
    # y ticks
    geom_segment(data = y.tick_frame, 
                 aes(x = zero - tick_sz, xend = zero + tick_sz, 
                     y = ticks, yend = ticks))+
    xlab( TeX('$SDE^*') ) +    scale_x_continuous(position = "top") +
    ylab(TeX('ME^*')) + 
    theme_classic() +
    labs(fill = NULL) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"))+
    theme(axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          text=element_text(size=16),#,  family="Palatino"),
          axis.ticks.x = element_blank(), axis.ticks.length.x = unit(-0.2, "cm"),
          axis.title.y = element_text(margin = margin(t = 0, r = -25, b = 0, l = 0), hjust = 0.02),
          legend.position= c(1,0.5),
          legend.justification="right",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-20), 
          legend.text.align = 0, 
          legend.title = element_text(vjust = 3)) +
    # labels
    geom_text(data=x.lab_frame, aes(x=lab, y=zero, label=lab),vjust=1.5, size = 4) +
    geom_text(data=y.lab_frame[-1,], aes(x=zero, y=lab, label=lab),hjust=1.5, size = 4)
  
  base.plot + 
    # labels
    geom_point(aes(color = colvar), shape = 19, size = 7)  +
    geom_point(aes(color = colvar), color = 'black', pch = 21, size = 7)  + 
    labs(color = colorval.name) +
    scale_color_viridis(option = "A") +
    {if(label == TRUE)geom_label_repel(aes(label = model),
                                       box.padding   = 0.35, 
                                       point.padding = 0.5,
                                       segment.color = 'grey50', 
                                       size = 5)}  +
    guides(colour = guide_colourbar(order = 1))
  
}

