DatasetSummary <- function(destination, pca){
  d.summary <- NULL
  filelist_path <- list.files(path=destination, pattern='.csv', full.names=T)
  filelist <- list.files(path=destination, pattern='.csv', full.names=F)
  
  for (i in 1:length(filelist_path)) {
    dataz <- read.csv(filelist_path[i], sep=';')
    dataz[,1] <- as.factor(dataz[,1])
    nindivs <- nlevels(dataz[,1])
    ncalls <- as.numeric(table(dataz[,1]))
    minncalls <- min(ncalls)
    maxncalls <- max(ncalls)
    nvars <- ncol(dataz[-1])
    nPCAs <- NA
    
    ### if pca=T, dataset is scaled and centered and variables are subjected to pca
    if(pca){
      # remove variables that have no variation
      eigs <- prcomp(dataz[,-1])$sdev^2
      varexplained <- eigs / sum(eigs)
      keepcols <- which(varexplained>0.0001)+1

      dataz <- calcPCAscaled(dataz)
      dataz <- data.frame(dataz[,1], dataz[,keepcols])
      nPCAs <- length(keepcols)
      plot(dataz[,2], dataz[,3], pch=c(0,1,2,3,4)[dataz[,1]], main=paste(filelist[i]), xlab='PC1', ylab='PC2', cex=1.5, cex.lab=2, cex.main=2)      } else {}

    HS <- calcHSnpergroup(dataz)[2]
    DS <- calcDS(dataz)
    HM <- calcHM(dataz)
    MI <- calcMI(dataz)
    datazpred <- data.frame(w=nindivs, x=round(mean(ncalls),0), y=DS)
    HSest <- predict(DStoHSloess, datazpred)
    datazpred <- data.frame(w=nindivs, x=round(mean(ncalls),0), z=HS)
    DSest <- predict(HStoDSloess, datazpred)
    
    d.summary <- rbind(d.summary, data.frame(filelist[i], nindivs, minncalls, 
                                             maxncalls, nvars, nPCAs, HS, DS, HM, MI, HSest, DSest))
  }
  
  row.names(d.summary) <- NULL
  return(d.summary)
}

