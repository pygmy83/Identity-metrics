ItterateDatasetsUni <- function(idRange, iRange, oRange, it, varmean) {
  result <- NULL
  for (id in 1:length(idRange)) {
    for (i in 1:length(iRange)) {
      ptm <- proc.time()  
      for (o in 1:length(oRange)) {
        for (it in itRange) {
          temp <- GenerateUnivariate(iRange[i], oRange[o], varmean, idRange[id]) 
          temp[,1] <- factor(temp[,1])
          grouping <- temp[,1]
          x <- temp[,2]
  
          Fvalue <- summary(aov(x~grouping))[[1]][['F value']][1]
          Pvalue <- summary(aov(x~grouping))[[1]][['Pr(>F)']][1]
          PICbwtotal <- calcPICbetweenmeans(temp)
          PICbwmeans <- calcPICbetweentot(temp)
          HSntot <- calcHSntot(temp)[2]
          HSngroups <- calcHSngroups(temp)[2]
          HSnpergroup <- calcHSnpergroup(temp)[2]
          HSvarcomp <- calcHSvarcomp(temp)
    
          result <- rbind(result, data.frame(iRange[i], oRange[o], idRange[id], it, Fvalue, Pvalue, 
                                             PICbwtotal, PICbwmeans,
                                             HSntot, HSngroups, HSnpergroup, HSvarcomp))
        }
      }
      
      print(i); print(max(iRange))
      print(proc.time()-ptm)
      
    }
  }
  
  names(result)[1:3] <- c('i', 'o', 'id')
  row.names(result) <- NULL
  return(result)
}
