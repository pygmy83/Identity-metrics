# test settings
iRange <- c(5,10,15) # value range for number of individuals:
oRange <- c(4,8,12) #, 12, 16, 20) # value range for number of observations = calls
pRange <- c(2,4) #, 6, 8, 10)  # value range for number of parameters = pcs
itRange <- c(1:5)# value range for number of iterations
idRange <- c(0.01, 1) # value range for individuality in data
covRange <- c(0,1) # value range for covariance
dataz <- ItterateDatasets(covRange,idRange,iRange,oRange,pRange,itRange) 

# settings for an initial analysis
iRange <- c(5, 10, 15, 20, 25, 30, 35, 40) # value range for number of individuals:
oRange <- c(4, 8, 12, 16, 20) #, 12, 16, 20) # value range for number of observations = calls
pRange <- c(2, 4, 6, 8, 10) #, 6, 8, 10)  # value range for number of parameters = pcs
itRange <- c(1:20)# value range for number of iterations
idRange <- c(0.01, 1, 2.5, 5, 10) # value range for individuality in data
covRange <- c(0, 0.25, 0.5, 0.75, 1) # value range for covariance 
dataz <- ItterateDatasets(covRange,idRange,iRange,oRange,pRange,itRange) 

# analysis to get HS values up to about 10 for DS and HS relationship
iRange <- c(4, 8, 12, 16, 20, 25, 30, 35, 40, 50, 75, 100) #, 15, 20, 25, 30, 35, 40, 45, 50) # range for number of individuals:
oRange <- c(4, 8, 12, 16, 20, 30, 40, 50, 75, 100) #, 12, 16, 20) # range for number of observations = calls
pRange <- c(2, 4, 6, 8, 10) #, 6, 8, 10)  # range for number of parameters = pcs
itRange <- c(1:10)# range for number of iterations
idRange <- c(0.1, 0.25, 0.5, 0.75, 1, 1.33, 1.66, 2)
covRange <- c(0)
dataz <- ItterateDatasets(covRange,idRange,iRange,oRange,pRange,itRange) 

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
            
            for (it in itRange) {
              #temp <- GenerateDataset(i, o, p, individualityRange[ir])
              #temp <- GenerateDataset2(iRange[i], oRange[o], pRange[p], covRange[cov], idRange[id])
              temp <- GenerateMultivariateDataset(iRange[i], oRange[o], pRange[p], covRange[cov], idRange[id])
              eigs <- prcomp(temp[,-1])$sdev^2
              varexplained <- eigs / sum(eigs)
              temp <- calcPCA(temp)
              #print(qplot(temp[,2], temp[,3], colour=temp[,1], main=paste0(cov, '-', it)))
              
              DS <- calcDS(temp)
              HS <- calcHSnpergroup(temp)[2]
              datazpred <- data.frame(w=iRange[i], x=oRange[o], y=DS)
              HSestLoess <- predict(m1, datazpred)
              HSestLin <- predict(m3, data.frame(y=DS))
              HSvarcomp <- calcHSvarcomp(temp)
              #HSvarcomp <- calcHSvarcomp(temp, varexplained)[1]
              MI <- calcMI(temp)
              #AMI <- MI[2]
              #NMI <- MI[1]
              HM <- calcHM(temp)
              #Fsel1 <- NA #calcFtype1(temp,2) # this is type I SS
              #Fsel2 <- calcFtype2(temp,2) # this is type II SS
              #ICC <- calcICC2(temp) # lmer method
              #auc <- calcAUC(temp)
              datazpred <- data.frame(w=iRange[i], x=oRange[o], z=HS)
              DSest <- predict(m1estDS, datazpred) # m1red - m1redestDS
              
              dataz <- rbind(dataz, data.frame(covRange[cov],idRange[id],iRange[i], oRange[o], pRange[p], it, 
                                                 DS, HS, HSestLoess, HSestLin, HSvarcomp, MI, HM, DSest))
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
  names(dataz)[1:5] <- c('cov', 'id', 'i', 'o', 'p')
  dataz <- round(dataz,3)
  return(dataz)
}


write.table(dataz, 'clipboard-100000', sep='\t', row.names=F)

# plot figures for Fig XX 'Multivariate identity metrics'
varstoplotX <- c(1:5)
varstoplotXnames <- names(dataz)[1:5]
varstoplotY <- c(7:10)
varstoplotYnames <- names(dataz)[7:10]
for (i in varstoplotX) {
  for (j in varstoplotY) {
    y <- dataz[,j]
    x <- dataz[,i]
    plotmeans(y~x, bars=T,n.label=F,barcol='black',ylab=varstoplotYnames[j], ylim=c(min(y), max(y)),xlab=varstoplotXnames[i], cex=2, cex.lab=2.5, cex.main=2.5, cex.axis=2)
    summary(aov(y~as.factor(x), data=dataz))
    TukeyHSD(aov(y~as.factor(x), data=dataz))
  }
}

y <- dataz$MI
x <- dataz$o
plot(y~x)


plotmeans(auc~i, data=nodata, bars=T,n.label=F, bar.col='black')
plot(temp[,2], temp[,3], col=temp$id, pch=20)
boxplot(temp[,2] ~ temp[,1])




#### comparison of real and estimated HS (loess) #####
HSest <- na.omit(nodata$HSestDSIO)
HSreal <- na.omit(nodata$HS)

HSest <- HSest[which(HSreal<8)]
HSreal <- HSreal[which(HSreal<8)]

summary(lm(HSest~HSreal))
plot(HSest~HSreal, xlim=c(0,15), ylim=c(0,15), ylab='HSest', col=rgb(0,100,0,10,maxColorValue=255), pch=16)
abline(0,1, col='grey', lwd=2, lty=2)
new.data <- predict(lm(HSest~HSreal),interval="prediction")  
S <- sqrt((sum((new.data[,1] - HSest)^2)) / length(HSest))
PI <- new.data[1,1] - new.data[1,2]
print(paste('Standard error of estimate:', round(S, 2)))
print(paste('Prediction interval:', round(PI, 2)))

lines(HSreal,new.data[,1],col="red",lwd=2)
lines(HSreal,new.data[,2],col="grey",lwd=2)
lines(HSreal,new.data[,3],col="grey",lwd=2)


#### comparison of real and estimated HS (linear) #####
HSest <- na.omit(nodata$HSestLin)
HSreal <- na.omit(nodata$HS)

HSest <- HSest[which(HSreal<8)]
HSreal <- HSreal[which(HSreal<8)]

summary(lm(HSest~HSreal))
plot(HSest~HSreal, xlim=c(0,10), ylim=c(0,10), ylab='HSest', col=rgb(0,100,0,10,maxColorValue=255), pch=16)
abline(0,1, col='grey', lwd=2, lty=2)
new.data <- predict(lm(HSest~HSreal),interval="prediction")  
S <- sqrt((sum((new.data[,1] - HSest)^2)) / length(HSest))
PI <- new.data[1,1] - new.data[1,2]
print(paste('Standard error of estimate:', round(S, 2)))
print(paste('Prediction interval:', round(PI, 2)))

lines(HSreal,new.data[,1],col="red",lwd=2)
lines(HSreal,new.data[,2],col="grey",lwd=2)
lines(HSreal,new.data[,3],col="grey",lwd=2)



#### comparison of real and estimated DS (loess) #####
DSest <- na.omit(nodata$DSest)
DSreal <- na.omit(nodata$DiscScore)

DSest <- na.omit(nodata$DSest[which(nodata$HS<8)])
DSreal <- na.omit(nodata$DiscScore[which(nodata$HS<8)])

summary(lm(DSest~DSreal))
plot(DSest~DSreal, xlim=c(0,1.2), ylim=c(0,1.2), ylab='DSest', col=rgb(0,100,0,10,maxColorValue=255), pch=16)
abline(0,1, col='grey', lwd=2, lty=2)
new.data <- predict(lm(DSest~DSreal),interval="prediction")  
S <- sqrt((sum((new.data[,1] - DSest)^2)) / length(DSest))
PI <- new.data[1,1] - new.data[1,2]
print(paste('Standard error of estimate:', round(S, 2)))
print(paste('Prediction interval:', round(PI, 2)))

lines(DSreal,new.data[,1],col="red",lwd=2)
lines(DSreal,new.data[,2],col="grey",lwd=2)
lines(DSreal,new.data[,3],col="grey",lwd=2)



#### comparison of real and estimated DS (linear) #####
#DSest <- na.omit(nodata$DSestLin)
#DSreal <- na.omit(nodata$DiscScore)

DSest <- na.omit(nodata$DSestLin[which(nodata$HS<8)])
DSreal <- na.omit(nodata$DiscScore[which(nodata$HS<8)])

summary(lm(DSest~DSreal))
plot(DSest~DSreal, xlim=c(0,1.2), ylim=c(0,1.2), ylab='DSest', col=rgb(0,100,0,10,maxColorValue=255), pch=16)
abline(0,1, col='grey', lwd=2, lty=2)
new.data <- predict(lm(DSest~DSreal),interval="prediction")  
S <- sqrt((sum((new.data[,1] - DSest)^2)) / length(DSest))
PI <- new.data[1,1] - new.data[1,2]
print(paste('Standard error of estimate:', round(S, 2)))
print(paste('Prediction interval:', round(PI, 2)))

lines(DSreal,new.data[,1],col="red",lwd=2)
lines(DSreal,new.data[,2],col="grey",lwd=2)
lines(DSreal,new.data[,3],col="grey",lwd=2)



###############################################

pairs(~DiscScore+i+o,data=nodata)
scatterplot.matrix(~DiscScore+i+o, data=nodata)

xdensity <- ggplot(nodata, aes(HS, fill=id)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00', '#999999', '#E69F00', '#999999')) + 
  theme(legend.position = "none")
xdensity


#library(rgl)
#library(car)
scatter3d(HS ~ i + o, data=nodata, fit="smooth")
dataz<-temp[temp$id==0.01,]
plot3d(dataz$i, dataz$o, dataz$HS)

individuality50ind.results <- nodata
individuality0.results <- nodata
individuality2.results <- nodata
individuality4.results <- nodata
individuality6.results <- nodata
individuality8.results <- nodata

plotmeans(DiscScore ~ o, bars=T, ylim=c(0, 1.2), xlab='number of calls', 
          ylab='discrimination score (10 indivs)', pch=0, data=individuality0.results)
plotmeans(DiscScore ~ o, bars=T, data=individuality2.results, pch=1, add=T)
plotmeans(DiscScore ~ o, bars=T, data=individuality4.results, pch=2, add=T)
plotmeans(DiscScore ~ o, bars=T, data=individuality6.results, pch=3, add=T)
plotmeans(DiscScore ~ o, bars=T, data=individuality8.results, pch=4, add=T)
legend(40,1.2, c("individuality=0", "individuality=2", "individuality=4", 'individuality=6', 'individuality=8'), pch=c(0,1,2,3,4))

plotmeans(HSall ~ o, bars=T, ylim=c(-1, 12), xlab='number of calls', 
          ylab='HS (10 indivs)', pch=0, data=individuality0.results)
plotmeans(HSall ~ o, bars=T, data=individuality2.results, pch=1, add=T)
plotmeans(HSall ~ o, bars=T, data=individuality4.results, pch=2, add=T)
plotmeans(HSall ~ o, bars=T, data=individuality6.results, pch=3, add=T)
plotmeans(HSall ~ o, bars=T, data=individuality8.results, pch=4, add=T)
legend(7,12, c("individuality=0", "individuality=2", "individuality=4", 'individuality=6', 'individuality=8'), pch=c(0,1,2,3,4))
