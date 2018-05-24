calcPICbetweentot <- function(temp) {
  temp[,1] <- as.factor(temp[,1])
  nvars <- length(names(temp))-1
  PIC <- rep(NA, nvars)
  
  for (i in 1:nvars){
    grouping <- temp[,1]
    x <- temp[,i+1]
    idmean <- tapply(x, grouping, mean)
    idsd <- tapply(x, grouping, sd)
    idn <- tapply(x, grouping, length)
    idwcv <- idsd / idmean
    #idwcv2 <- 100 * (1 + 1/(4*idn)) * (idsd / idmean)  # correction for low sample sizes
    meanWCV <- mean(idwcv, na.rm=T)
    #meanWCV2 <- mean(idwcv2, na.rm=T)   # correction for low sample size
    
    BCV <- sd(x) / mean(x) ### ### Coefficient of variation between individuals calculated from all samples
    
    PIC[i] <- BCV / meanWCV
    #PIC[i] <- (100 * BCV * (1 + 1/(4*length(x)))) / meanWCV2 # correction for low sample size
    
  }
  
  names(PIC) <- names(temp[-1])
  return(PIC)
}
