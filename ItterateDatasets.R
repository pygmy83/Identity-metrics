ItterateDatasets <- function(covRange,idRange,iRange,oRange,pRange,itRange){
  dataz <- NULL
  tottime <- proc.time()
  
  for (cov in 1:length(covRange)){
    
    for (id in 1:length(idRange)){
      
      for (i in 1:length(iRange)) {
        #i <- 7
        ptm <- proc.time()  
        
        for (o in 1:length(oRange)) {
          #o <- 10  
          
          for (p in 1:length(pRange)) {
            #p <- 4  
            
            for (it in 1:length(itRange)) {
              #temp <- GenerateDataset(i, o, p, individualityRange[ir])
              #temp <- GenerateDataset2(iRange[i], oRange[o], pRange[p], covRange[cov], idRange[id])
              temp <- GenerateMultivariateDataset(iRange[i], oRange[o], pRange[p], covRange[cov], idRange[id])
              eigs <- prcomp(temp[,-1])$sdev^2
              varexplained <- eigs / sum(eigs)
              temp <- calcPCA(temp)
              #print(qplot(temp[,2], temp[,3], colour=temp[,1], main=paste0(cov, '-', it)))
              
              DS <- calcDS(temp)
              HS <- calcHSnpergroup(temp)[2]
              MI <- calcMI(temp)
              HM <- calcHM(temp)

              datazpred <- data.frame(w=iRange[i], x=oRange[o], y=DS)
              HSestLoess <- predict(DStoHSloess, datazpred)
              HSestLin <- predict(DStoHSlin, data.frame(y=DS))
              
              datazpred <- data.frame(w=iRange[i], x=oRange[o], z=HS)
              DSestLoess <- predict(HStoDSloess, datazpred) 
              DSestLin <- predict(HStoDSlin, data.frame(z=HS))
              
              dataz <- rbind(dataz, data.frame(covRange[cov],idRange[id],iRange[i], oRange[o], pRange[p], itRange[it], 
                                                 DS, HS, MI, HM, HSestLoess, HSestLin, DSestLoess, DSestLin))
              #datasetname <- paste(i,o,p,it)
              #write.table(temp, file=paste0('dataset', datasetname,'.csv'), sep=';', row.names = F)
            } #it
            
          } #p
        } #o
        
        #write.csv(dataz, file='simulationbackup.csv', row.names = F)
        print(paste('cov:',cov, 'of', length(covRange))); print(paste('id:',id, 'of', length(idRange))); print(paste('individuals:',i, 'of', length(iRange)))
        print(proc.time()-ptm)
        
      } #i
    } #ir
  }# cov
  print('total time:'); print(proc.time()-tottime)
  row.names(dataz) <- NULL
  names(dataz)[1:6] <- c('cov', 'id', 'i', 'o', 'p', 'it')
  dataz <- round(dataz,3)
  return(dataz)
}

#write.table(dataz, 'clipboard-100000', sep='\t', row.names=F)

#y <- dataz$HM
#x <- dataz$i
#plotmeans(y~x, bars=T,n.label=F,barcol='black',
#          yaxt ='n', xaxt='n',
#          cex=2, cex.lab=2, cex.axis=2)

