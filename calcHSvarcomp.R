calcHSvarcomp <- function (temp, sumHS=T) {
  
    nvars <- ncol(temp)
    vars <- names(temp)
    #n <- nrow(temp)
    #indivs <- levels(as.factor(temp[,1]))
    #nindiv <- length(indivs)
    #npergroup <- nrow(temp) / nindiv
    #fvalues <- rep(NA, nvars)
    # Pr <- rep(NA, nvars) P-value is not calculated for mixed model
    HS <- rep(NA, nvars) 

    randomVar <- rep(NA,nvars)
    residVar <- rep(NA,nvars)
    totalVar <- rep(NA,nvars)
    
    for (k in 2:nvars) {
      modelFormula <- paste0(vars[k], '~1+(1|', vars[1], ')')
      lmer.m1 <- lmer(as.formula(modelFormula), data=temp)
      randomVar[k] <- (as.data.frame(VarCorr(lmer.m1))[1,4])
      residVar[k] <- as.data.frame(VarCorr(lmer.m1))[2,4]
      totalVar[k] <- randomVar[k]+residVar[k]
      HS[k] <- log2((totalVar[k] / residVar[k]))
    }
    
    # Pr <- round(Pr, 3)
    HS <- round(HS, 2)
    result <- data.frame(vars,HS)

    if (sumHS==T) {
      result <- c(sum(result$HS,na.rm=T))
      names(result) <- c('HS for all vars')
      return (result)
    } else return (result)
}
