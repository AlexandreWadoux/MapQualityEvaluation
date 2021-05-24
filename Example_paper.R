########################################
########################################
#simulate field
# Gaussian process with mean zero and exponential covariance function that has sill of 5 and a range of 10
# add a linear spatial trend using the coordinates
# add a multiplicative term 

#define discretisation grid: AW: for the paper I used 120, not 70
grid <- expand.grid(x1 = seq(0, 70, length.out = 70),
                    x2 = seq(0, 70, length.out = 70))

#compute spatial trend; x1 and x2 ares used as covariate, the linear model has an intercept
grid$mu <- 5 +  0.1*grid$x1 + 0.05*grid$x2

#define covariance function for simulation of residuals
covfun <- function(sill, range, Dist) {
  sill * exp(-Dist / range)
}

#compute matrix with distances between simulation nodes
distx<-outer(grid$x1,grid$x1,FUN="-")
disty<-outer(grid$x2,grid$x2,FUN="-")
Dist<-sqrt(distx^2+disty^2)
#Dist <- as.matrix(dist(grid[c('x1', 'x2')]^2))

#compute matrix with mean covariances
sill <- 5
range <- 10
C <- covfun(sill, range, Dist = Dist)

#simulate values for residuals by Cholesky decomposition
set.seed(31415)
Upper <- chol(C)
G <- rnorm(n=nrow(grid),0,1) #simulate random numbers from standard normal distribution
grid$residuals <- crossprod(Upper,G)

# add the trend and the Cholesky decomposed values
grid$sim <- grid$mu+grid$residuals

# add the multiplicative effect
grid$y <- grid$sim*0.3

#add residuals to trend
grid <- grid[,c('x1', 'x2', 'y')]

library(ggplot2)
ggplot(grid) + geom_tile(aes(x1, x2, fill = y)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

mean(grid$y)
var(grid$y)

# the gridL is 0-70, we will keep 0-50. AW: for the paper the grid is from 0-120, and we keep keep 0-100
gridL <- grid
grid$x1[grid$x1>50] <- NA
grid$x2[grid$x2>50] <- NA
grid <- na.omit(grid)

########################################
########################################
# 1- positional error
gridL$x1[gridL$x1<20] <- NA
gridL$x2[gridL$x2<20] <- NA
gridL <- na.omit(gridL)
library(ggplot2)
ggplot(gridL) + geom_tile(aes(x1, x2, fill = y)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

grid$y_posError <- gridL$y

########################################
########################################
# 1- reverse simulated field

grid$y_rev <- rev(grid$y)
ggplot(grid) + geom_tile(aes(x1, x2, fill = y_rev)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

########################################
########################################
# 2- average simulated field

grid$y_mean <- mean(grid$y)
ggplot(grid) + geom_tile(aes(x1, x2, fill = y_mean)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()


########################################
########################################
# 3- smoothed simulated field

# we use a moving window average
## convert to raster
library(sp)
coordinates(grid) <- ~x1+x2
gridded(grid) <- TRUE
library(raster)
r <- stack(grid)
plot(r)

# make the smoothing
r$y_smoothed <- focal(r$y, w = matrix(1, nc=5, nr=5), fun=mean, na.rm=T, pad = T)
plot(r$y_smoothed)

# convert raster back to data.frame
grid <- as.data.frame(r, xy = TRUE)
colnames(grid)[1:2] <-  c('x1', 'x2')

ggplot(grid) + geom_tile(aes(x1, x2, fill = y_smoothed)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

########################################
########################################
# 4- biased simulated field

grid$y_bias <- grid$y + 0.5 # AW: in the paper I used 1

ggplot(grid) + geom_tile(aes(x1, x2, fill = y_bias)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

########################################
########################################
# 5- RF predicted field

# make a RF model with a single tree (it is not equivalent to a regression tree)
library(ranger)
model <- ranger(y ~ x1 + x2 , data = grid, num.trees = 1)
grid$y_RFpred <- predict(object = model, as.data.frame(grid))$predictions

ggplot(grid) + geom_tile(aes(x1, x2, fill = y_RFpred)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

library(DescTools)
CCC(grid$y, grid$y_RFpred)$rho.c
cor(grid$y, grid$y_RFpred)^2

# make a regression tree model 
library(rpart)
model <- rpart(y ~ x1 + x2 , data = grid)
grid$y_RTpred <- predict(object = model, grid)

ggplot(grid) + geom_tile(aes(x1, x2, fill = y_RTpred)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

library(DescTools)
CCC(grid$y, grid$y_RTpred)$rho.c
cor(grid$y, grid$y_RTpred)^2


########################################
########################################
# 6- only high or low values from the mean
# take the mean of the simulated field
mean_ySim <- mean(grid$y)

# TEST 1
# take the values higher than the 75% quantile, otherwise return the mean
grid$y_high <- NA
for(i in 1:nrow(grid)){
  if(grid$y[i] < quantile(grid$y, c(.75)) ){grid$y_high[i] <- mean_ySim}else{grid$y_high[i] <- grid$y[i]}
}

ggplot(grid) + geom_tile(aes(x1, x2, fill = y_high)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

library(DescTools)
CCC(grid$y, grid$y_high)$rho.c
cor(grid$y, grid$y_high)^2

# TEST 2
# take the values lower than the 25% quantile, otherwise return the mean
grid$y_low <- NA
for(i in 1:nrow(grid)){
  if(grid$y[i] > quantile(grid$y, c(.25))){grid$y_low[i] <- mean_ySim}else{grid$y_low[i] <- grid$y[i]}
}

ggplot(grid) + geom_tile(aes(x1, x2, fill = y_low)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()


ras2 <- grid 
coordinates(ras2) <- ~x1+x2
gridded(ras2) <- TRUE
ras2 <- stack(ras2)
names(ras2)[1] <- 'y_original'


library(rasterVis)
library(viridis)
library(RColorBrewer)
miat = seq(0, 9, length.out = 100)
myTheme <- rasterTheme(region = viridis(11),
                       strip.background = list(col = 'transparent'), 
                       axis.line = list(col = "transparent"),
                       layout.heights=list(xlab.key.padding = 1)) 
pdf(file= './Simulated_case_maps.pdf',family="Palatino", width = 13, height = 6)             
levelplot(ras2,
          at = miat,
          names.attr = c('Reference', 'Positional shift', 'Reversed', 'Mean', 'Smoothed', 'Negative bias', 
                         'Random forest', 'Regression tree', 'Upper quartile', 'Lower quartile' ),
          layout=c(5, 2),
          scales=list(x=list(rot=30, cex=0.8), y=list(rot=30, cex=0.5), col = 'transparent'),
          colorkey=list(space="top"),
          par.settings = myTheme
) 
dev.off()


########################################
########################################
# the sun diagram
library(tdr)
library(ggrepel)
library(magrittr)
library(latex2exp)
library(viridis)
source('./R_sun_function.R')
source('./map_quality_indices.R')
models <- list(y_rev = grid$y_rev,
               y_posError = grid$y_posError,
               y_mean = grid$y_mean,
               y_smoothed = grid$y_smoothed, 
               y_bias = grid$y_bias,
               y_RFpred = grid$y_RFpred,
               y_RTpred = grid$y_RTpred, 
               y_high = grid$y_high, 
               y_low = grid$y_low)
names(models) <- c('Reversed', 'Positional error', 'Mean', 'Smoothed', 'Negative bias', 
                   'Random forest', 'Regression tree', 'Upper quartile', 'Lower quartile' )
obser <- grid$y

# additional values to be plotted as colour on the diagram, the r
model_MEC <- plyr::laply(models, .fun = function(x,y){eval(x,y)$NSE}, y = grid$y)
colorval <- model_MEC

pdf(file= './Simulated_case_solar.pdf',family="Palatino", width = 12, height = 11)             
gg_sun(mods = models, 
       obs = obser,
       colorval = colorval,
       colorval.name = 'MEC',
       label = TRUE, 
       x.axis_begin = -1.1,
       x.axis_end = 1.1,
       y.axis_end = 1.9,
       by = 0.1)
dev.off()


########################################
########################################
# the target diagram
library(tdr)
library(ggrepel)
library(magrittr)
library(latex2exp)
library(viridis)
source('./R_target_function.R')

#AW: to be checked: the sign must be 1 or -1, otherise sd(obs)-sd(pred) has sign 0 for y_rev and y_rev is in the center because it has no bias
models <- list(y_rev = grid$y_rev,
               y_posError = grid$y_posError,
               y_mean = grid$y_mean,
               y_smoothed = grid$y_smoothed, 
               y_bias = grid$y_bias,
               y_RFpred = grid$y_RFpred,
               y_RTpred = grid$y_RTpred, 
               y_high = grid$y_high, 
               y_low = grid$y_low)
names(models) <- c('Reversed', 'Positional error', 'Mean', 'Smoothed', 'Negative bias', 
                   'Random forest', 'Regression tree', 'Upper quartile', 'Lower quartile' )
obser <- grid$y

# additional values to be plotted as colour on the diagram, the r
model_cors <- plyr::laply(models, .fun = function(x,y){cor(x,y)}, y = grid$y)
# account for possible sd=0 (mean field), then cor = 0
model_cors[is.na(model_cors)] <- 0
colorval <- model_cors

pdf(file= './Simulated_case_target.pdf',family="Palatino", width = 13, height = 12)             
gg_target(mods = models, 
          obs = obser,
          colorval = colorval,
          colorval.name = 'Correlation',
          label = TRUE, 
          axis_begin  =-2.2,
          axis_end   = 2.2,
          by = 0.1)
dev.off()

# Test: Eq. 14 Jolliff
M07 <- sqrt(1+(0.7^2) -2*0.7^2)
M1 <- sqrt(1+(1^2) -2*1^2)


########################################
########################################
# the taylor diagram
library(ggforce)
library(ggrepel)
library(latex2exp)
source('R_taylor_function2.R')

models <- list(y_rev = grid$y_rev,
               y_posError = grid$y_posError,
               y_mean = grid$y_mean,
               y_smoothed = grid$y_smoothed, 
               y_bias = grid$y_bias,
               y_RFpred = grid$y_RFpred,
               y_RTpred = grid$y_RTpred, 
               y_high = grid$y_high, 
               y_low = grid$y_low)
names(models) <- c('Reversed', 'Positional error', 'Mean', 'Smoothed', 'Negative bias', 
                   'Random forest', 'Regression tree', 'Upper quartile', 'Lower quartile' )
obser <- grid$y

pdf(file= './Simulated_case_taylor.pdf',family="Palatino", width = 12, height = 7)             
gg_taylor(mods = models, 
          obs = obser, 
          label = TRUE)
dev.off()


##########################
##########################
# compute map quality indices
source('map_quality_indices.R')
tab.res <- plyr::laply(models, .fun = eval, y = obser)
rownames(tab.res) <-  c('Reversed', 'Positional error', 'Mean', 'Smoothed', 'Negative bias', 
                        'Random forest', 'Regression tree', 'Upper quartile', 'Lower quartile' )
tab.res

##########################
##########################
# test: compute distance from reference point, see Eq 8 in Jollif et al., (2009)
# Calculate the statistics
model_cors <- plyr::laply(models, .fun = cor, y = obser)
# if cor = NA because of sd = 0 (when using the mean of the data as model), change to cor = 0, this is realistic, see discussion is Rstats
model_cors[is.na(model_cors)] <- 0
model_stds <- plyr::laply(models, .fun = sd)
# normalize
obs_std <- 1
model_stds <- model_stds/sd(obser, na.rm = T)

data.frame(model = c('y_rev','y_posi', 'y_mean', 'y_smoothed', 'y_bias', 'y_RFpred', 'y_RTpred', 'y_high', 'y_low'), 
           distRefCenter = sqrt(1+model_stds^2 - 2*model_stds*model_cors))



















