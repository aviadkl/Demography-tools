## create Sprague multiplier data
## Insert data with five-year age groups with any last open group (75+, 80+, 85+ etc.)

sprague <- function(x,last.age){
  
  breaks=c(seq(0,last.age,5),150)

  Ns <- x

  group_num <- length(Ns)

  last_pop_plus <- Ns[group_num]

  multipliers <- data.frame(G1 = c(0.3616, 0.264, 0.184, 0.12, 0.0704,
            0.0336, 0.008, -0.008, -0.016, -0.0176,
            -0.0128,  -0.0016, 0.0064, 0.0064, 0.0016,
            rep(0,10)),
    G2 = c(-0.2768,  -0.096, 0.04, 0.136, 0.1968,
            0.2272, 0.232, 0.216, 0.184, 0.1408,
            0.0848, 0.0144, -0.0336, -0.0416, -0.024,
            -0.0144,  -0.008, 0, 0.008, 0.0144,
            0.0176, 0.016, 0.008, -0.008, -0.0336),
    G3 = c(0.1488, 0.04, -0.032, -0.072, -0.0848,
           -0.0752, -0.048, -0.008, 0.04, 0.0912,
            0.1504, 0.2224, 0.2544, 0.2224, 0.1504,
            0.0912, 0.04, -0.008, -0.048,  -0.0752,
            -0.0848, -0.072,  -0.032, 0.04, 0.1488),
    G4 = c(-0.0336, -0.008, 0.008,  0.016, 0.0176,
            0.0144,  0.008, 0, -0.008, -0.0144,
            -0.024,  -0.0416, -0.0336,  0.0144, 0.0848,
            0.1408,  0.184,  0.216, 0.232, 0.2272,
            0.1968, 0.136,  0.04, -0.096, -0.2768),
    G5 = c(0,  0, 0,  0, 0,
            0,  0, 0, 0, 0,
            0.0016,  0.0064, 0.0064, -0.0016, -0.0128,
            -0.0176, -0.016, -0.008, 0.008, 0.0336,
            0.0704, 0.12,  0.184, 0.264, 0.3616)
  )
  multipliers <- cbind(groups=rep(c("lowest","low","normal","high", "highest"),
                                  each=5), multipliers)

  # 27 years:
  ## correspond to n3 since it is
  ## in group 25,26,27,28,29
  infoGroup <- function(n, mult=multipliers, mybreaks=breaks, popN=Ns){
    ## from a five years group, which one of the five years:
    groups <- (n) %% 5
    ## extreme group or normal:
    if(n < 5){
      tab <- subset(mult, subset=groups=="lowest")
    } else if(n >= 5 & n < 10){
      tab <- subset(mult, subset=groups=="low")
    } else if(n >= (last.age-5) & n < last.age){
      tab <- subset(mult, subset=groups=="high")
    } else if(n >= (last.age-1)){
      tab <- subset(mult, subset=groups=="highest")
    } else{
      tab <- subset(mult, subset=groups=="normal")
    }
    ## my groups:
    ng <- cut(n, mybreaks, right=FALSE)
    mygroup <- which(levels(ng) %in% ng)
    ## cohort
    rowsm <- tab[groups+1,2:ncol(tab)]
    if(mygroup==1){
      s <- seq(mygroup,mygroup+4,1)
    } else if(mygroup==2){
      s <- seq(mygroup-1, mygroup+3,1)
    } else if(mygroup==group_num){
      s <- seq(mygroup-4,mygroup,1)
    } else if(mygroup==group_num-1){
      s <- seq(mygroup-3,mygroup+1,1)
    } else{
      s <- seq(mygroup-2,mygroup+2, 1)
    }


    cohort <- sum(rowsm * popN[s])
    return(cohort)
  }

   cohorts <- sapply(0:(last.age - 1), infoGroup)

  ## last age group+
   cohorts <- c(cohorts, last_pop_plus)
   plus <- paste(as.character(last.age), "+")
   names(cohorts) <- c(0:(last.age - 1),plus)
   return(cohorts)
}

	