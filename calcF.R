#### type I sum of squares for all variables in a dataset #####
calcFtype1 <- function(temp) {
  temp[,1] <- as.factor(temp[,1])
  nvars <- length(names(temp))-1
  f <- rep(NA, nvars)
  
  for (i in 1:nvars){
    modelFormula <- paste(names(temp)[i+1], '~', names(temp)[1])
    f[i] <- summary(aov(as.formula(modelFormula), data=temp))[[1]][["F value"]][[1]]
  }
  names(f) <- names(temp[-1])
  return(f)
}

#### type 2 sum of squares as originally in Beecher 1989, for all variables in a dataset ####
#### needs 'car' library ####
library(car)
calcFtype2 <- function(temp) {
  temp[,1] <- as.factor(temp[,1])
  nvars <- length(names(temp))-1
  f <- rep(NA, nvars)
  
  for (i in 1:nvars){
    modelFormula <- paste(names(temp)[i+1], '~', names(temp)[1])
    f[i] <- Anova(lm(as.formula(modelFormula), data=temp), type='II')[['F value']][1]
  }
  names(f) <- names(temp[-1])
  return(f)
}




