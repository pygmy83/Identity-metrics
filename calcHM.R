### 'calcHM2' uses three other functions 'calcMeanVec', 'calcDistT', and 'calcDistW' (included belllow)

calcHM <- function(temp){
  temp[,1] <- as.factor(temp[,1])
  Y <- temp[,-1]
  MeanVec <- calcMeanVec(Y)
  distT <- calcDistT(temp)
  distW <- calcDistW(temp)

  g <- length(levels(temp[,1]))
  npergroup <- length(temp[,1]) / g
  FM <- ((npergroup-1) / (g-1)) * ((distT - g*(distW)) / distW)
  HM <- log2(sqrt(((FM+npergroup-1) / npergroup)))
  
  return(HM)
}


### calcMeanVec - calculates centroid of measurments for each individual
calcMeanVec <- function(Y){
  nvars <- ncol(Y)
  MeanVec <- rep(NA, nvars)
  
  for (i in 1:nvars) {
    MeanVec[i] <- mean(Y[,i]) 
  }
  return(MeanVec)
}


### calcDistT - calculates sum of distances between each sample and centroid of all samples
calcDistT <- function(temp){
  Y <- temp[-1]
  MeanVec <- calcMeanVec(Y)
  MeanVecDist <- rep(NA, nrow(Y))
  for (i in 1:nrow(Y)){
    MeanVecDist[i] <- sum((MeanVec-Y[i,])^2) * (1/length(MeanVec))
  }
  TotalDist <- sum(MeanVecDist)
  return(TotalDist)
}


### calcDistW - calculates sum of distances between each group centroid and samples within the group and returns the average within group distance
calcDistW <- function (temp) {
  indivs <- levels(temp[,1])
  nindivs <- length(indivs)
  DistW <- rep(NA, nindivs)
    for (i in 1:nindivs){
      temp.id <- temp[temp[,1]==indivs[i],]
      DistW[i] <- calcDistT(temp.id)
    }
  return(mean(DistW))
}


