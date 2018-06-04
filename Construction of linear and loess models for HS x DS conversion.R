# build the model
# the data used to build model was generated through ItterateDatasets function with following settings:
# analysis to get HS values up to about 10 for DS and HS conversions
# iRange <- c(5, 15, 20, 25, 30, 35, 40) # range for number of individuals:
# oRange <- c(4, 8, 12, 16, 20) # range for number of observations = calls
# pRange <- c(2, 4, 6, 8, 10) #, 6, 8, 10)  # range for number of parameters = pcs
# itRange <- c(1:10)# range for number of iterations
# idRange <- c(0.1, 0.25, 0.5, 0.75, 1, 1.33, 1.66, 2)
# covRange <- c(0)
# dataz <- ItterateDatasets(covRange,idRange,iRange,oRange,pRange,itRange) 

dataz <- read.delim('data-HSxDSconversion-model building data.csv', sep=';', header=T) # load data for model building

w <- dataz$i
x <- dataz$o
y <- dataz$DS
z <- dataz$HS
  
m1 <- loess(z ~ w + x + y, span=.7, degree=2) # DS to HS loess
m2 <- loess(y ~ w + x + z, span=.7, degree=2) # HS to DS loess
m3 <- lm(z ~ y) # DS to HS linear
m4 <- lm(y ~ z) # HS to DS linear
  



# test the model
# separate dataset was used to test the model; the dataset was build again with the same settings as the previous dataset

dataz <- read.delim('data-HSxDSconversion-model testing data.csv', header=T, sep=';') 


pdf("mygraph.pdf", width=7, height=7)
par(mfrow=c(2,2))



#### comparison of real and estimated DS (linear) #####
DSest <- na.omit(dataz$DSestLin)
DSreal <- na.omit(dataz$DS)
#DSest <- na.omit(dataz$DSestLin[which(dataz$HS<8)])
#DSreal <- na.omit(dataz$DS[which(dataz$HS<8)])

summary(lm(DSest~DSreal))
plot(DSest~DSreal, xlim=c(0,1.15), xlab='DS',ylim=c(0,1.15), ylab='DS (linear estimate)', col=rgb(0,0,0,10,maxColorValue=255), pch=16)
abline(0,1, col='grey', lwd=1, lty=2)
new.data <- predict(lm(DSest~DSreal),interval="prediction")  
S <- sqrt((sum((new.data[,1] - DSest)^2)) / length(DSest))
PI <- new.data[1,1] - new.data[1,2]
print(paste('Standard error of estimate:', round(S, 2)))
print(paste('Prediction interval:', round(PI, 2)))

lines(DSreal,new.data[,1],col="black",lwd=1)
lines(DSreal,new.data[,2],col="grey",lwd=1)
lines(DSreal,new.data[,3],col="grey",lwd=1)



#### comparison of real and estimated DS (loess) #####
DSest <- na.omit(dataz$DSestLoess)
DSreal <- na.omit(dataz$DS)
#DSest <- na.omit(dataz$DSestLoess[which(dataz$HS<8)])
#DSreal <- na.omit(dataz$DS[which(dataz$HS<8)])

summary(lm(DSest~DSreal))
plot(DSest~DSreal, xlim=c(0,1.15), xlab='DS', ylim=c(0,1.15), ylab='DS (loess estimate)', col=rgb(0,0,0,10,maxColorValue=255), pch=16)
abline(0,1, col='grey', lwd=1, lty=2)
new.data <- predict(lm(DSest~DSreal),interval="prediction")  
S <- sqrt((sum((new.data[,1] - DSest)^2)) / length(DSest))
PI <- new.data[1,1] - new.data[1,2]
print(paste('Standard error of estimate:', round(S, 2)))
print(paste('Prediction interval:', round(PI, 2)))

lines(DSreal,new.data[,1],col="black",lwd=1)
lines(DSreal,new.data[,2],col="grey",lwd=1)
lines(DSreal,new.data[,3],col="grey",lwd=1)



#### comparison of real and estimated HS (linear) #####
HSest <- na.omit(dataz$HSestLin)
HSreal <- na.omit(dataz$HS)
#HSest <- HSest[which(HSreal<8)]
#HSreal <- HSreal[which(HSreal<8)]

summary(lm(HSest~HSreal))
plot(HSest~HSreal, xlim=c(0,15), xlab='HS', ylim=c(0,15), ylab='HS (linear estimate)', col=rgb(0,0,0,10,maxColorValue=255), pch=16)
abline(0,1, col='grey', lwd=1, lty=2)
new.data <- predict(lm(HSest~HSreal),interval="prediction")  
S <- sqrt((sum((new.data[,1] - HSest)^2)) / length(HSest))
PI <- new.data[1,1] - new.data[1,2]
print(paste('Standard error of estimate:', round(S, 2)))
print(paste('Prediction interval:', round(PI, 2)))

lines(HSreal,new.data[,1],col="black",lwd=1)
lines(HSreal,new.data[,2],col="grey",lwd=1)
lines(HSreal,new.data[,3],col="grey",lwd=1)



#### comparison of real and estimated HS (loess) #####
HSest <- na.omit(dataz$HsestLoess)
HSreal <- na.omit(dataz$HS)
#HSest <- HSest[which(HSreal<8)]
#HSreal <- HSreal[which(HSreal<8)]

summary(lm(HSest~HSreal))
plot(HSest~HSreal, xlim=c(0,15), xlab='HS', ylim=c(0,15), ylab='HS (loess estimate)', col=rgb(0,0,0,10,maxColorValue=255), pch=16)
abline(0,1, col='grey', lwd=1, lty=2)
new.data <- predict(lm(HSest~HSreal),interval="prediction")  
S <- sqrt((sum((new.data[,1] - HSest)^2)) / length(HSest))
PI <- new.data[1,1] - new.data[1,2]
print(paste('Standard error of estimate:', round(S, 2)))
print(paste('Prediction interval:', round(PI, 2)))

lines(HSreal,new.data[,1],col="black",lwd=1)
lines(HSreal,new.data[,2],col="grey",lwd=1)
lines(HSreal,new.data[,3],col="grey",lwd=1)

dev.off()

  





