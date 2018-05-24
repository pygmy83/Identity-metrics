calcHSnpergroup <- function (temp, sumHS=T){
  
  nvars <- ncol(temp)
  vars <- names(temp)
  n <- nrow(temp)
  indivs <- levels(as.factor(temp[,1]))
  nindiv <- length(indivs)
  npergroup <- nrow(temp) / nindiv
  fvalues <- rep(NA, nvars)
  Pr <- rep(NA, nvars)
  HS <- rep(NA, nvars) 
  
  for (k in 2:nvars) {
    modelFormula <- paste(vars[k], '~', vars[1])
    fvalues [k] <- summary(aov(as.formula(modelFormula), data=temp))[[1]][["F value"]][[1]]
    Pr[k] <- summary(aov(as.formula(modelFormula), data=temp))[[1]][["Pr(>F)"]][[1]]
    #nindiv[k] <- nlevels(temp[,1])
    HS [k] <- log2(sqrt((fvalues[k]+(npergroup-1)) / npergroup)) #npergroup
  }
  Pr <- round(Pr, 3)
  HS <- round(HS, 2)
  result <- data.frame(vars,Pr,HS)
  
  if (sumHS==T) {
    result <- c(sum(result$HS[result$Pr<0.05],na.rm=T), sum(result$HS,na.rm=T))
    names(result) <- c('HS for significant vars', 'HS for all vars')
    return (result)
  } else return (result)
}
