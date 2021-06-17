eval <- function(x, y){

    # mean error
    ME <- round(mean(y - x, na.rm = TRUE), digits = 2)
    
    # root mean square error
    RMSE <-   round(sqrt(mean((y - x)^2, na.rm = TRUE)), digits = 2)
    
    # root mean absolute error
    MAE <-   round(mean(abs(y - x), na.rm = TRUE), digits = 2)
    
    # Pearson's correlation squared
    r2 <-  round((cor(x, y, method = 'pearson', use = 'pairwise.complete.obs')^2), digits = 2)
    
    # MEC
    SSE <- sum((y - x) ^ 2, na.rm = T)
    SST <- sum((y - mean(y, na.rm = T)) ^ 2, na.rm = T)
    NSE <- round((1 - SSE/SST), digits = 2)
    
    # concordance correlation coefficient
    n <- length(x)
    sdx <- sd(x, na.rm = T)
    sdy <- sd(y, na.rm = T)
    r <- stats::cor(x, y, method = 'pearson', use = 'pairwise.complete.obs')
    # scale shift
    v <- sdx / sdy
    sx2 <- var(x, na.rm = T) * (n - 1) / n
    sy2 <- var(y, na.rm = T) * (n - 1) / n
    # location shift relative to scale
    u <- (mean(x, na.rm = T) - mean(y, na.rm = T)) / ((sx2 * sy2)^0.25)
    Cb <- ((v + 1 / v + u^2)/2)^-1
    rCb <- r * Cb
    rhoC <- round(rCb, digits = 2)
    
    Cb <- round(Cb, digits = 2)
    r <- round(r, digits = 2)
    
    # return the results
    evalRes <- data.frame(ME = ME, MAE = MAE, RMSE = RMSE, r = r, r2 = r2, NSE = NSE, rhoC = rhoC, Cb = Cb)
  
  return(evalRes)
}
