##### Life Table function - complete or abridged #####

# two functions: 1.life table by mx input.  2.life table by qx input. 

### qx = vector of survival probabilities
### Mx = vector of mortality rates
### age.interval = 1 or 5
### max.age = last age group, e.g 80+, 85+ etc.

# column 1 = nx
# column 2 = ax
# column 3 = Mx
# column 4 = qx
# column 5 = px
# column 6 = lx
# column 7 = dx
# column 8 = Lx
# column 9 = Tx
# column 10 = ex
# column 11 = Sx


life.table.mx <- function(mx, age.interval, max.age){
  
  if (age.interval==1){
    table <- matrix(data = NA, nrow = max.age+1, ncol = 11)
    table[,1] <- 1
    
  } else{
    table <- matrix(data = NA, nrow = (max.age/5 + 2), ncol = 11)
    table[,1] <- c(1,4,rep(5,max.age/5))
  }
  
  table[,2] <- c(0.1,rep(0.5,nrow(table)-1))
  table[,3] <- mx
  table[,4] <- table[,1]*table[,3]/(1+table[,1]*(1-table[,2])*table[,3])
  table[nrow(table),4] <- 1
    
  table[,5] <- 1 - table[,4]
  
  table[1,6] <- 100000
  for (i in 2:nrow(table)){
    table[i,6] <- table[i-1,5]*table[i-1,6]
  }
  
  
  table[,7] <- table[,6]*table[,4]
  
  
  for (i in 1:(nrow(table)-1)){
    table[i,8] <- table[i,1]*(table[i+1,6]+(table[i,2]*table[i,7]))
  }
  table[nrow(table),8] <- table[nrow(table),6]/table[nrow(table),3] 
   
  
  for (i in 1:nrow(table)){
    table[i,9] <- sum(table[i:nrow(table),8])
  }
  
  table[,10] <- table[,9]/table[,6]
  
  if (age.interval==1){
    for (i in 1:(nrow(table)-1)){
      table[i,11] <- table[i+1,8]/table[i,8]
      if (i==nrow(table)-1) table[i,11] <- table[i+1,9]/table[i,9]
    }
  }else{
    table[2,11] <- table[3,8]/sum(table[1:2,8])#survival 0-4
    table[nrow(table)-1,11] <- table[nrow(table),9]/table[nrow(table)-1,9]
    for (i in 3:(nrow(table)-2)){
      table[i,11] <- table[i+1,8]/table[i,8]
    }
  }

  return(table)
  
}



life.table.qx <- function(qx, age.interval, max.age, last.ax=0.5){
  
  if (age.interval==1){
    table <- matrix(data = NA, nrow = max.age+1, ncol = 11)
    table[,1] <- 1
    
  } else{
    table <- matrix(data = NA, nrow = (max.age/5 + 2), ncol = 11)
    table[,1] <- c(1,4,rep(5,max.age/5))
  }
  
  table[,2] <- c(0.1,rep(0.5,nrow(table)-1))
  table[nrow(table),2] <- last.ax;
  
  table[,4] <- qx
  table[nrow(table),4] <- 1
  
  table[,3] <- table[,4]/(table[,1]-(table[,1]-table[,2])*table[,4])
  
  table[,5] <- 1 - table[,4]
  
  table[1,6] <- 100000
  for (i in 2:nrow(table)){
    table[i,6] <- table[i-1,5]*table[i-1,6]
  }
  
  
  table[,7] <- table[,6]*table[,4]
  
  
  for (i in 1:(nrow(table)-1)){
    table[i,8] <- table[i,1]*(table[i+1,6]+(table[i,2]*table[i,7]))
  }
  table[nrow(table),8] <- table[nrow(table),2]/table[nrow(table),6] 
  
  
  for (i in 1:nrow(table)){
    table[i,9] <- sum(table[i:nrow(table),8])
  }
  
  table[,10] <- table[,9]/table[,6]
  table[nrow(table),10] <- last.ax
  table[nrow(table),3] <- 1/(table[nrow(table),10]) #last value of mx, that can now be completed
  
  if (age.interval==1){
    for (i in 1:(nrow(table)-1)){
      table[i,11] <- table[i+1,8]/table[i,8]
      if (i==nrow(table)-1) table[i,11] <- table[i+1,9]/table[i,9]
    }
  }else{
    table[2,11] <- table[3,8]/sum(table[1:2,8])#survival 0-4
    table[nrow(table)-1,11] <- table[nrow(table),9]/table[nrow(table)-1,9]
    for (i in 3:(nrow(table)-2)){
      table[i,11] <- table[i+1,8]/table[i,8]
    }
  }
  
  return(table)
  
}