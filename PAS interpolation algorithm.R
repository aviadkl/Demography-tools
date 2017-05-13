### algorithm that interpolates mortality rates from a range of two Mx vectors ###
### source of the algorithm is the PAS spreadsheets by the US census bureau ###

### currently works only for single age data ###

# e0.min, e0.max - life expectancy at age 0 for both vectors
# eLastAge.min, eLastAge.max - life expectancy at last age for both vectors
# mx.min, mx.max - mortality rates vectors for both vectors
# e0.target - life expectancy at age 0 for which Mx's are required



Mx.interpolate <- function(e0.min, e0.max, eLastAge.min, eLastAge.max, mx.min, mx.max, e0.target){
  
  source("Q:\\Aviad\\R code\\life table function.R")
  
  mx.min <- as.matrix(mx.min)
  mx.max <- as.matrix(mx.max)
  
  mx.min.lifetable <- life.table.mx(mx=mx.min, age.interval=1, max.age=(NROW(mx.min)-1))
  qx1 <- as.matrix(mx.min.lifetable[,4])
  
  mx.max.lifetable <- life.table.mx(mx=mx.max, age.interval=1, max.age=(NROW(mx.max)-1))
  qx2 <- as.matrix(mx.max.lifetable[,4])
  
  ln.qx1 <- log(qx1)
  ln.qx2 <- log(qx2)
  
  # Estimated Life Tables function
  est.life <- function(factor){
    
    factor.inv <- 1 - factor
    
    est.life.table <- matrix(data = NA, nrow = NROW(qx1), ncol = 4)
    
    est.life.table[,1] <- c(0.1,rep(0.5,NROW(qx1)-1))   #nx*ax
    est.life.table[NROW(qx1),1] <- factor*eLastAge.min + factor.inv*eLastAge.max   #last age weight
    
    est.life.table[,2] <- exp(factor*ln.qx1 + factor.inv*ln.qx2)   #qx
    est.life.table[NROW(qx1),2] <- 1 #last qx
    
    est.life.table[1,3] <- 100000 #lx
    for (i in 2:NROW(qx1)){
      est.life.table[i,3] <- (1-est.life.table[i-1,2])*est.life.table[i-1,3]
    }
    
    for (i in 1:(NROW(qx1)-1)){  #Lx
      est.life.table[i,4] <- est.life.table[i+1,3] + est.life.table[i,1]*(est.life.table[i,3]-est.life.table[i+1,3])
    }
    est.life.table[NROW(qx1),4] <- est.life.table[NROW(qx1),1]*est.life.table[NROW(qx1),3]  #last Lx
    
    T0 <- sum(est.life.table[,4])
    e0.est <- T0/100000
    
    # quadratic constants
    a2 <- (factor.inv*(e0.max-e0.est)+factor*(e0.min-e0.est))/((e0.est*e0.est-e0.max*e0.max)*(e0.min-e0.est)+(e0.min*e0.min-e0.est*e0.est)*(e0.max-e0.est))
    a1 <- (a2*(e0.est*e0.est-e0.max*e0.max)-factor)/(e0.max-e0.est)
    a0 <- -a1*e0.max-a2*e0.max*e0.max
    
    output <- list(e0.est, a2, a1, a0, est.life.table[,2])
    return(output)
    
  }
  
  
  
  # First Estimated Life Table
  factor1 <- (e0.max-e0.target)/(e0.max-e0.min)
  est1 <- est.life(factor1)
  
  # Second Estimated Life Table
  est2.a2 <- est1[[2]]
  est2.a1 <- est1[[3]]
  est2.a0 <- est1[[4]]
  factor2 <- est2.a0+est2.a1*e0.target+est2.a2*e0.target*e0.target
  est2 <- est.life(factor2) 
   
  # Third Estimated Life Table
  est3.a2 <- est2[[2]]
  est3.a1 <- est2[[3]]
  est3.a0 <- est2[[4]]
  factor3 <- est3.a0+est3.a1*e0.target+est3.a2*e0.target*e0.target
  est3 <- est.life(factor3) 
  
  # Last Estimated Life Table
  est4.a2 <- est3[[2]]
  est4.a1 <- est3[[3]]
  est4.a0 <- est3[[4]]
  factor4 <- est4.a0+est4.a1*e0.target+est4.a2*e0.target*e0.target
  est4 <- est.life(factor4) 
  
  final.qx <- as.matrix(est4[[5]])
  last.ax <- factor4*eLastAge.min+(1-factor4)*eLastAge.max     #important for last mx calculation
  final.lifetable <- life.table.qx(qx=final.qx, age.interval=1, max.age=(NROW(mx.min)-1), last.ax=last.ax)
  final.mx <- as.matrix(final.lifetable[,3])
  
  return(final.mx)
  
}