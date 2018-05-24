destination <- './empirical datasets examples'
temp <- DatasetSummary(destination, pca=T)
write.table(temp, paste0(destination, '/datasetsummary.txt'), sep='\t')

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
      # remove constant variables
      eigs <- prcomp(dataz[,-1])$sdev^2
      varexplained <- eigs / sum(eigs)
      keepcols <- which(varexplained>0.0001)+1
      #keepcols <- c(2,3)
      
      dataz <- calcPCAscaled(dataz)
      #print(filelist[i]); print(names(dataz)[keepcols])
      dataz <- data.frame(dataz[,1], dataz[,keepcols])
      nPCAs <- length(keepcols)
      #print(qplot(dataz[,2], dataz[,3], colour=dataz[,1], main=paste(filelist[i])))
      plot(dataz[,2], dataz[,3], pch=c(0,1,2,3,4)[dataz[,1]], main=paste(filelist[i]), xlab='PC1', ylab='PC2', cex=1.5, cex.lab=2, cex.main=2)      } else {}
      #plot(dataz[,2], dataz[,3], col=dataz[,1], pch=20, main=paste(filelist[i]), xlab='PC1', ylab='PC2')      } else {}
  
    HS <- calcHSnpergroup(dataz)[2]
    DS <- calcDS(dataz)
#    HM <- calcHM2(dataz)
#    MI <- calcMI(dataz)
#    AMI <- MI[2]
#    NMI <- MI[1]
#    MI <- MI[3]
#    ICC <- calcICC2(dataz)
#    AUC <- as.numeric(calcAUC(dataz))
    datazpred <- data.frame(w=nindivs, x=round(mean(ncalls),0), y=DS)
    #HSestDSIO <- predict(m1, datazpred)
#    HSestDSIO <- predict(m1, datazpred)
    #HSestDSIO <- predict(m1exp, datazpred)
#    NMIest <- calcNMIest(dataz[,1], DS)
    datazpred <- data.frame(w=nindivs, x=round(mean(ncalls),0), z=HS)
#    DSest <- predict(m1estDS, datazpred) # m1red - m1redestDS
    
    d.summary <- rbind(d.summary, data.frame(filelist[i], nindivs, minncalls, 
                                             maxncalls, nvars, nPCAs, HS, DS))
  }
  
  #d.summary$HSrank <- rank(datasetsummary$HS, ties.method= "average")
  #d.summary$DSrank <- round(rank(datasetsummary$DS, ties.method= "average"), 0)
  #d.summary$HMrank <- rank(datasetsummary$HM, ties.method= "average")
  #d.summary$AUCrank <- rank(datasetsummary$AUC, ties.method= "average")
  
  #return(d.summary[order(datasetsummary$HSrank, decreasing=T),])
  return(d.summary)
}

