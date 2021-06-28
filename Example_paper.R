
############
############  load libraries
############

library(ggplot2)
library(ranger)
library(rpart)
library(DescTools)
library(sp)
library(raster)
library(rasterVis)
library(viridis)
library(RColorBrewer)
library(tdr)
library(ggrepel)
library(magrittr)
library(latex2exp)
library(ggforce)

############
############  Initialize script
############

root_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_dir)

############
############ Simulate fields 
############
# Reference field 
## grid of size 50*50 cells
## we simulate a realisation of a map composed of a linear spatial trend superimposed on a Gaussian random field
## trend has an intercept of 5 and slope parameter of 0.1 for the x-axis and 0.05 for the y-axis
## the Gaussian random field has a mean of zero and a covariance given by C(h) = 5 exp(-h/10), where h is the lag distance
## this map is multiplied by a factor of 0.3 to obtain the ultimate reference map

# Nine modifications of the reference map are created to give nine maps with predictions: 
## 1- a map positionally shifted by 20 units in the x-direction and the y-direction,
## 2- a map where the x- and y-coordinates are reversed, 
## 3- a map containing the mean value of the reference map at each coordinate, 
## 4- a smoothed map obtained by a moving average window of size 5*5 units,  
## 5- a negatively biased map obtained by adding a value of 1 to the reference map, 
## 6- a map predicted by a random forest model with a single tree and default parameters from the R package ranger using the reference map as calibration data and the x- and y-coordinates as predictors,
## 7- a map predicted by a single regression tree with default parameters from the R package rpart using the reference map as calibration data and the x- and y-coordinates as predictors, 
## 8- a map of the upper quartile, made by assigning the mean of the reference map to the values lower than the upper quartile and value of the reference map otherwise,
## 9- a map of the lower quartile, made by assigning the mean of the reference map to the values higher than the lower quartile and the value of the reference map otherwise.

########################################
########################################
# Reference field 
#define discretisation grid: AW: for the paper I used 120, not 70 (we remove 20 after, so 50 or 100)
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

# plot reference map
ggplot(grid) + geom_tile(aes(x1, x2, fill = y)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

# mean and variance of the reference map
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
grid$y_posError <- gridL$y
ggplot(gridL) + geom_tile(aes(x1, x2, fill = y)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

########################################
########################################
# 2- reverse simulated field

grid$y_rev <- rev(grid$y)
ggplot(grid) + geom_tile(aes(x1, x2, fill = y_rev)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

########################################
########################################
# 3- average simulated field

grid$y_mean <- mean(grid$y)
ggplot(grid) + geom_tile(aes(x1, x2, fill = y_mean)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

########################################
########################################
# 4- smoothed simulated field

# we use a moving window average
## convert to raster
coordinates(grid) <- ~x1+x2
gridded(grid) <- TRUE
r <- stack(grid)
plot(r)

# make the smoothing
r$y_smoothed <- focal(r$y, w = matrix(1, nc = 5, nr = 5), fun = mean, na.rm = T, pad = T)
plot(r$y_smoothed)

# convert raster back to data.frame
grid <- as.data.frame(r, xy = TRUE)
colnames(grid)[1:2] <-  c('x1', 'x2')

ggplot(grid) + geom_tile(aes(x1, x2, fill = y_smoothed)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

########################################
########################################
# 5- biased simulated field
grid$y_bias <- grid$y + 0.5 # AW: in the paper I used 1

ggplot(grid) + geom_tile(aes(x1, x2, fill = y_bias)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

########################################
########################################
# 6- RF predicted field
# make a RF model with a single tree (it is not equivalent to a regression tree)
model <- ranger(y ~ x1 + x2 , data = grid, num.trees = 1)
grid$y_RFpred <- predict(object = model, as.data.frame(grid))$predictions

ggplot(grid) + geom_tile(aes(x1, x2, fill = y_RFpred)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

########################################
########################################
# 7- CART predicted field
# make a regression tree model 
model <- rpart(y ~ x1 + x2 , data = grid)
grid$y_RTpred <- predict(object = model, grid)

ggplot(grid) + geom_tile(aes(x1, x2, fill = y_RTpred)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

########################################
########################################
# 8- upper quartile map
# take the mean of the simulated field
mean_ySim <- mean(grid$y)

# take the values higher than the 75% quantile, otherwise return the mean
grid$y_high <- NA
for(i in 1:nrow(grid)){
  if(grid$y[i] < quantile(grid$y, c(.75)) ){grid$y_high[i] <- mean_ySim}else{grid$y_high[i] <- grid$y[i]}
}

ggplot(grid) + geom_tile(aes(x1, x2, fill = y_high)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

########################################
########################################
# 9- lower quartile map
# take the values lower than the 25% quantile, otherwise return the mean
grid$y_low <- NA
for(i in 1:nrow(grid)){
  if(grid$y[i] > quantile(grid$y, c(.25))){grid$y_low[i] <- mean_ySim}else{grid$y_low[i] <- grid$y[i]}
}

ggplot(grid) + geom_tile(aes(x1, x2, fill = y_low)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

########################################
########################################
# convert the dataframe to a raster 
ras2 <- grid 
coordinates(ras2) <- ~x1+x2
gridded(ras2) <- TRUE
ras2 <- stack(ras2)
names(ras2)[1] <- 'y_original'

########################################
########################################
# plot
miat = seq(0, 6, length.out = 100)
myTheme <- rasterTheme(region = viridis(11),
                       strip.background = list(col = 'transparent'), 
                       axis.line = list(col = "transparent"),
                       layout.heights=list(xlab.key.padding = 1)) 

jpeg(file= './Simulated_case_maps.jpg',family="Palatino", width = 26, height = 12, units = "cm", res = 1200)            
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


############
############  solar diagram
############

# load the functions
source('./R_solar_function.R')
source('./map_quality_indices.R')

# make a list of vectors (map values)
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

# vector of true values or observations
obser <- grid$y

# additional values to be plotted as colour on the diagram, the MEC (NSE)
model_MEC <- plyr::laply(models, .fun = function(x,y){eval(x,y)$NSE}, y = grid$y)
colorval <- model_MEC
 
jpeg(file= './Simulated_case_solar.jpg', family="Palatino", width = 25, height = 21, units = "cm", res = 2500)             
gg_solar(mods = models, 
         obs = obser,
         colorval = colorval,
         colorval.name = 'MEC',
         label = TRUE, 
         x.axis_begin = -1.1,
         x.axis_end = 1.1,
         y.axis_end = 1.9,
         by = 0.1)
dev.off()


############
############  target diagram (jolliff et al. 2009)
############

# load the function
source('./R_target_function.R')

# make a list of vectors (map values)
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

# vector of true values or observations
obser <- grid$y

# additional values to be plotted as colour on the diagram, the MEC
model_MEC <- plyr::laply(models, .fun = function(x,y){eval(x,y)$NSE}, y = grid$y)
colorval <- model_MEC

jpeg(file= './Simulated_case_target.jpg', family="Palatino", width = 24, height = 22, units = "cm", res = 2500)             
gg_target(mods = models, 
          obs = obser,
          colorval = colorval,
          colorval.name = 'MEC',
          label = TRUE, 
          axis_begin  =-2.2,
          axis_end   = 2.2,
          by = 0.1)
dev.off()

# Test: Eq. 14 Jolliff
M07 <- sqrt(1+(0.7^2) -2*0.7^2)
M1 <- sqrt(1+(1^2) -2*1^2)


############
############  Taylor diagram (Taylor, 2001)
############

# load the function
source('R_taylor_function.R')

# make a list of vectors (map values)
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

# vector of true values or observations
obser <- grid$y

jpeg(file= './Simulated_case_taylor.jpg', family="Palatino", width = 24, height = 14, units = "cm", res = 2500)            
gg_taylor(mods = models, 
          obs = obser, 
          label = TRUE)
dev.off()


############
############  Map quality indices
############

# load the function
source('map_quality_indices.R')

# create a table and compute the indices
tab.res <- plyr::laply(models, .fun = eval, y = obser)
rownames(tab.res) <-  c('Reversed', 'Positional error', 'Mean', 'Smoothed', 'Negative bias', 
                        'Random forest', 'Regression tree', 'Upper quartile', 'Lower quartile' )
# view the table
tab.res

############
############  Distance between reference and other points in the target diagram
############

# compute distance from reference point, see Eq 8 in Jolliff et al. (2009)

# Calculate the statistics
model_cors <- plyr::laply(models, .fun = cor, y = obser)

# if cor = NA because of sd = 0 (when using the mean of the data as model), change to cor = 0, this is realistic, see discussion is Rstats
model_cors[is.na(model_cors)] <- 0
model_stds <- plyr::laply(models, .fun = sd)

# standardize
obs_std <- 1
model_stds <- model_stds/sd(obser, na.rm = T)

# store the results in a data frame
data.frame(model = c('y_rev','y_posi', 'y_mean', 'y_smoothed', 'y_bias', 'y_RFpred', 'y_RTpred', 'y_high', 'y_low'), 
           distRefCenter = sqrt(1+model_stds^2 - 2*model_stds*model_cors))



















