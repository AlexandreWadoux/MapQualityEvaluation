gg_taylor <- function(mods, obs, label){
  library(ggforce)
  library(ggrepel)
  
  # Calculate the statistics
  model_cors <- plyr::laply(mods, .fun = cor, y = obs)
  # if cor = NA because of sd = 0 (when using the mean of the data as model), change to cor = 0, this is realistic, see discussion is Rstats
  model_cors[is.na(model_cors)] <- 0
  model_stds <- plyr::laply(mods, .fun = sd)
  obs_std <- 1
  model_stds <- model_stds/sd(obs, na.rm = T)
  xlabTitle <-  TeX('Standardized standard deviation $\\sigma^* $')
  
  model_points_df <- data.frame(Model = names(mods), Cor = model_cors, Std = model_stds)
  
  std_max <- ceiling((max(c(obs_std, model_stds, 2))))
  std_major <- seq(0, std_max, length.out = 5)
  
  semicircle_df <- plyr::adply(std_major, 1, .fun = function(x){data.frame(x = x * cos(seq(0,pi,pi/1000)),
                                                                           y = x * sin(seq(0,pi,pi/1000)),
                                                                           label = x)})
  cor_major <- c(seq(-1, 1, 0.1),-0.95, 0.95, -0.99, 0.99)
  lines_df <- plyr::adply(cor_major, 1, .fun = function(x){data.frame(xend = 1*max(std_major) * cos(acos(x)),
                                                                      yend = 1*max(std_major) * sin(acos(x)),
                                                                      label = x)})
  
  # add the distance from the reference 
  x_0 = 1
  max_R = 2
  ref_r_min=0.25
  ref_r_max=2
  ref_contours = 8
  ref_r = seq(ref_r_min,ref_r_max,length.out = ref_contours)
  
  x_star = array(0,dim=ref_contours)
  x_star <- (ref_r^2 - max_R^2 - x_0^2)/(-2*x_0) 
  x_star[ref_r < (max_R-x_0) ] = x_0+ref_r[ref_r < (max_R-x_0) ] 
  x_star[x_star<0]=NaN  
  
  # Now determine the y points where we intersect
  y_star=sqrt(ref_r^2-(x_star-x_0)^2)
  y_star[ref_r < (max_R-x_0)]=0 
  
  angle_end=pi/2-atan(y_star/((x_star-x_0))) 
  angle_end[x_star < x_0] <- -pi/2+atan(y_star[x_star < x_0]/(x_0-x_star[x_star < x_0]))
  
  
  angle_begin <- -pi/2
  angle_begin[is.nan(angle_begin)] = -pi/2
  
  small_arcs <- data.frame(angle_begin,angle_end,ref_r)#
  
  # standard deviation around the reference point
  maxray <- 1.5 * std_max
  discrete <- seq(180, 0, by = -1)
  dat.circle <- data.frame(xcircle = numeric(), ycircle = numeric(), label = numeric())
  seq.sd <- seq.int(0, 2 * maxray, by = (maxray/12))[-1]
  
  # create the lines
  for (i in seq.sd) {
    xcircle <- 1 + (cos(discrete * pi/180) * i)
    ycircle <- sin(discrete * pi/180) * i
    labelc <- rep(i, length(xcircle))
    dat.circle <- rbind(dat.circle, cbind(xcircle, ycircle, labelc))
  }
  F <- approxfun(x=semicircle_df[semicircle_df$label == max(semicircle_df$label),]$x, 
                 semicircle_df[semicircle_df$label == max(semicircle_df$label),]$y)
  
  # remove the lines outside the plot
  for(i in 1:nrow(dat.circle)){
    if(dat.circle$ycircle[i] > F(dat.circle$xcircle[i]) | is.na(F(dat.circle$xcircle[i]))){
      dat.circle$ycircle[i] <- NA
      dat.circle$xcircle[i] <- NA}
  }    
  
  # prepare labels
  dat.circle.labelc <- data.frame(xcircle = numeric(), ycircle = numeric(), label = numeric())
  for (i in unique(dat.circle$labelc)){
    dat.sub <- dat.circle[dat.circle$labelc == i,]
    dat.circle.labelc <- rbind(dat.circle.labelc, dat.sub[10,])}
  
  base_plot <- 
    ggplot() +
    coord_equal()+
    ggthemes::theme_base() +
    geom_line(data = semicircle_df, aes(x = x, y = y, group = label), linetype = "solid", colour = "black", size = 0.6) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)))+
    scale_x_continuous(labels = abs)+
    geom_segment(data = lines_df, aes(x = 0, y = 0, xend = xend, yend = yend), linetype = "dashed", colour = "black") +
    geom_segment(data = lines_df, aes(x = min(xend), y = 0, xend = max(xend), yend = 0), colour = "black") +
    geom_line(data = dat.circle, aes(x = xcircle, y = ycircle, group = labelc), linetype = "dashed", colour = "red3", size = 0.6) +
    geom_label(data = dat.circle.labelc, aes(x = xcircle, y = ycircle, label = labelc),label.size = NA , fill = "white", vjust = 0, size = 5, colour = "red3") +
    geom_text(data = lines_df, aes(x = 1.07 * xend, y = 1.035 * yend, label = label), vjust = 0, size = 5, colour = "black") +
    theme_taylor(base_size = 12)+
    xlab(xlabTitle)+
    annotate("text", x=0, y=std_max+0.25, label= TeX("Correlation \\textit{r}"), size = 5)
  
  base_plot + 
    geom_point(data = model_points_df, aes(x = Std * cos(acos(Cor)), y = Std * sin(acos(Cor))), size = 6) + 
    geom_point(aes(x =1, y = 0), size = 3, colour ='red3') +
    {if(label == TRUE) geom_label_repel(data = model_points_df, aes(x = Std * cos(acos(Cor)), y = Std * sin(acos(Cor)), label = Model),
                                        box.padding   = 0.35, 
                                        point.padding = 0.8,
                                        segment.color = 'grey50', 
                                        size = 5)} + 
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          text=element_text(size=16,  family="Palatino"),
          panel.border = element_blank(),
          axis.line.x = element_blank(),
          panel.grid = element_blank())
}

theme_taylor <- function(base_size = 11){
  theme_bw(base_size = base_size) %+replace%
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          text=element_text(size=16,  family = "Palatino"),
          panel.border = element_blank(),
          axis.line.x = element_blank(),
          panel.grid = element_blank())
}
