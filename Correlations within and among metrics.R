
### prune different number of individuals and calls in each dataset ###################################
#######################################################################################################

destination <- './data simul'
outputdestination <- './data pruned'
callrange <- c(4:20)
indivrange <- c(5:40)
seed <- 1
DatasetPrune(destination, outputdestination, callrange, indivrange, seed)


DatasetPrune <- function(destination, outputdestination, callrange, indivrange, seed) {
  
  filelist_path <- list.files(path=destination, pattern='.csv', full.names=T)
  filelist <- list.files(path=destination, pattern='.csv', full.names=F)
  nfiles <- length(filelist)
  ncalls <- sample(callrange, nfiles,replace=T)
  nindivs <- sample(indivrange, nfiles,replace=T)
  
  
  for (i in 1:length(filelist_path)) {
    dataz <- read.csv(filelist_path[i], sep=';')
    dataz[,1] <- as.factor(dataz[,1])
    
    ### select individuals with more calls per individuals than 'ncalls'
    a <- names(which(table(dataz[,1])>=ncalls[i]))
    temp <- dataz[dataz[,1] %in% a, ]
    temp[,1] <- factor(temp[,1])
    
    ### check if there is more individuals than 'nindivs' left after
    if (length(a) < nindivs) stop(paste('there is less than', nindivs, 'individuals with more than', ncalls, 
                                        'calls for', filelist[i], 'dataset' ))
    
    ### take random sample of 'nindivs' individuals 
    ### take random sample of 'ncalls' from each individual
    set.seed(seed)
    indivs <- sample(levels(as.factor(temp[,1])), nindivs[i])
    rndid <- with(temp, ave(temp[,2], temp[,1], FUN=function(x) {sample.int(length(x))}))
    obs <- which(rndid <= ncalls[i])
    temp2 <- dataz[obs, ]
    temp3 <- temp2[temp2[,1] %in% indivs, ]
    temp3[,1] <- factor(temp3[,1])
    
    write.table(temp3, paste0(outputdestination,'/', filelist[i]), sep=';', row.names=F)
  }
}


### subsample repeatedly and calculate metrics ############################################################
###########################################################################################################

### summary of full datasets
destination <- './simulated datasets/'
summary.all <- DatasetSummary(destination,pca=T)
write.table(summary.all, paste0(destination, '/summaryall.txt'), sep='\t',row.names=F)

#temp <- summary.all[,7:12]
#write.table(cor(temp, method='spearman'), 'clipboard', sep='\t', row.names=T)


### Prune sequantially the full dataset and write summary
result <- NULL
iterations <- c(1:2)
tottime <- proc.time()
callrange <- c(2:10)
indivrange <- c(2:10)
seed <- 1
destination <- './simulated datasets/'
outputdestination <- './pruned datasets/'

for (i in iterations){
  ptm <- proc.time()
  do.call(file.remove, list(list.files(outputdestination, full.names = TRUE)))
  DatasetPrune(destination, outputdestination, callrange, indivrange, seed)
  
  destination <- './pruned datasets/' 
  temp <- DatasetSummary(destination,pca=T)
  temp$it <- i
  print("Number of individuals in each resampled datasets:")
  print(temp$nindivs)
  result <- rbind(result, temp)
  print(i)
  print(proc.time()-ptm)
}; print('total time:'); print(proc.time()-tottime)




### compute correlations between metrics from full dataset in 'summary.all' and partial datasets in 'result'
### - within metric consistency ############################################################################
############################################################################################################

result <- read.delim('data-consistency-result-simulated.csv', header=T, sep=';')
summary.all <- read.delim('data-consistency-summary.all-simulated.csv', header=T, sep=';')
result$it <- as.factor(result$it)
cormet <- 'pearson' # correlation method
a <- calcWithinConsist(summary.all, result, metrics, cormet)

calcWithinConsist <- function(summary.all, result, metrics, cormet){
  result.cor.within <- NULL
  for (i in 1:nlevels(result$it)){
    temp <- result[result$it==i,]  
    HScor <- cor(summary.all[,7], temp[,7], method=cormet, use="complete.obs")
    DScor <- cor(summary.all[,8], temp[,8], method=cormet, use="complete.obs")
    HMcor <- cor(summary.all[,9], temp[,9], method=cormet, use="complete.obs")
    MIcor <- cor(summary.all[,10], temp[,10], method=cormet, use="complete.obs")
    temp <- data.frame(HScor, DScor, HMcor, MIcor)
    result.cor.within <- rbind(result.cor.within, temp)
  } 
  return(result.cor.within)
}

boxplot(a, use.cols = T, ylim=c(0.75,1))
friedman.test(as.matrix(a))
posthoc.friedman.nemenyi.test(as.matrix(a)) # library(PMCMR)

t.test(result.cor.within[,2], result.cor.within[,8], paired=T)

##############################################################################################################



### compute correlations in 'result' - among metric correlations ################
result <- read.delim('data-consistency-result-empirical.csv', header=T, sep=';')
summary.all <- read.delim('data-consistency-summary.all-empirical.csv', header=T, sep=';')
result.cor.among <- NULL
result$it <- as.factor(result$it)

for (i in 1:nlevels(result$it)){
  temp <- result[result$it==i,]  
  DScor.all <- cor(summary.all[,8], temp[,8], method='pearson', use="complete.obs")
  DSDSestcor.all <- cor(summary.all[,8], temp[,12], method='pearson', use="complete.obs")
  DSDSestcor <- cor(temp[,8], temp[,12], method='pearson', use="complete.obs")
  HSHScor.all <- cor(summary.all[,7], temp[,7], method='pearson', use="complete.obs")
  HSHSestcor.all <- cor(summary.all[,7], temp[,11], method='pearson', use="complete.obs")
  HSHSestcor <- cor(temp[,7], temp[,11], method='pearson', use="complete.obs")

  temp <- data.frame(DScor.all, DSDSestcor.all, DSDSestcor, HSHScor.all, HSHSestcor.all, HSHSestcor)
  result.cor.among <- rbind(result.cor.among, temp)
}

summary(result.cor.among)

# DS and DSest correlations
boxplot(result.cor.among[,1], result.cor.among[,2], result.cor.among[,3], names=names(result.cor.among)[c(1,2, 3)])
a <- result.cor.among[,1:3]
friedman.test(as.matrix(a))
posthoc.friedman.nemenyi.test(as.matrix(a)) # library(PMCMR)


# HS and HSest correlations
boxplot(result.cor.among[,4], result.cor.among[,5], result.cor.among[,6], names=names(result.cor.among)[c(4,5,6)])
a <- result.cor.among[,4:6]
friedman.test(as.matrix(a))
posthoc.friedman.nemenyi.test(as.matrix(a)) # library(PMCMR)

