# replicate Searby Jouventin results
it <- 3 # number of iterations
nindivs <- 50 # number of individuals
ncalls <- 6 # number of calls per individual
nvars <- 25 # number of variables measured
popul <- 20 # number of populations
BI <- c(0.00, 0.07, 0.13, 0.19, 0.24, 0.29, 0.34, 0.38, 0.43, 0.46, 0.50, 0.53, 0.56, 0.59, 0.62, 0.65, 0.67, 0.69, 0.71, 0.73) # range of values for between individual variation
WI <- c(1.00,0.93,0.87,0.81,0.76,0.71,0.66,0.62,0.57,0.54,0.50,0.47,0.44,0.41,0.38,0.35,0.33,0.31,0.29,0.27) # range of values for within individual variation
temp <- GenerateDataset4(it, nindivs, ncalls, nvars, popul, BI, WI) # see function GenerateDatset4 below

HM <- c()
DS <- c()
HS <- c()
SDtot <- c()
SDwithin <- c()
SDbetween <- c()

# This loop calculates HS, DS and HM from the simulated dataset for each of the 20 populations and 3 iterations
for (i in 1:it){
  for (j in 1:popul) {
    a <- temp[temp$it==i,]
    a <- a[a$populid==j,]
    #write.table(a, 'clipboard-100000', sep='\t', row.names=F)
    a <- a[,5:30]
    HM <- c(HM, calcHM(a)) # vector with HM values
    DS <- c(DS, calcDS(a)) # vector with DS values
    HS <- c(HS, calcHSnpergroup(a)[2]) # vector with HS values
    SDtot <- c(SDtot, sd(a[,2]))
    SDwithin <- c(SDwithin, calcSDwithin (a[,1:2]))
    SDbetween <- c(SDbetween, calcSDbetween (a[,1:2]))
  }
}

HM2 <- 2^HM
HS2 <- 2^HS
dataz <- data.frame(HM, HM2,DS, HS, HS2)

#### This is to replicate relationship between DS and 2^HM as reported in Searby and Jouventin (2004)
a <- 1
x <- seq(1,2, by=0.01)
x0 <- 1.22
b <- -11.5
y <- a / (1+ (x/x0)^b)

#### Plot the results and compare with the relationship reported in Searby and Jouventin (2004, Fig 2b)
plot(DS~HM2, axes=F, xlim=c(0.8,2.4), ylim=c(-0.2,1.2), pch=16)
axis(side=1, at=seq(0.8, 2.4, by=0.2))
axis(side=2, at=seq(0, 1.2, by=0.2))
box()
points(x, y, pch=3, col='red') # add points illustrating results reported in Searby and Jouventin (2008)




### HS and HM in our simulated datasets with changes in covariance and number of variables
dataz <- read.delim('data-IterateDatasetsMulti.csv', header=T, sep=';')
plot(dataz$HS~dataz$HM, main='all data', xlab='HM', ylab='HS')

dataz <- dataz[dataz$cov==0,]
dataz <- dataz[dataz$i==40,]
dataz10 <- dataz[dataz$p==10,]
plot(dataz10$HS~dataz10$HM, main='covariance=0, number of individuals=40', xlab='HM', ylab='HS')
summary(lm(dataz10$HS~dataz10$HM))
abline(lm(dataz10$HS~dataz10$HM))

dataz4 <- dataz[dataz$p==4,]
points(dataz4$HM, dataz4$HS, col='grey')
summary(lm(dataz4$HS~dataz4$HM))
abline(lm(dataz4$HS~dataz4$HM), col='grey')



####################################################################
#### Function generating dataset as in Searby and Jouventin (2004)##
####################################################################

GenerateDataset4 <- function(it, nindivs, ncalls, nvars, popul, BI, WI){
  indcodes <- seq(from=1, to=nindivs, by=1)
  parammat <- matrix(rep(NA, it*popul*ncalls*nindivs*nvars), nrow=it*popul*nindivs*ncalls, ncol=nvars, byrow=F)
  
  for (p in 1:nvars) {
    paramvec <- c()
    id <- c()
    itnumber <- c()
    populid <- c()
    populWI <- c()
    populBI <- c()
    
    for (k in 1:it) {
      for (j in 1:popul) {
        indmeans <- rnorm(nindivs, 0, sqrt(BI[j]))
        for (i in 1:nindivs) {
          paramvec <- c(paramvec, rnorm(ncalls, indmeans[i], sqrt(WI[j])))
          id <- c(id, rep(indcodes[i], ncalls))
          populid <- c(populid, rep(j,ncalls))
          populWI <- c(populWI, rep(WI[j], ncalls))
          populBI <- c(populBI, rep(BI[j], ncalls))
          itnumber <- c(itnumber, rep(k,ncalls))
        } # indivs
      } # popul
    } #it
    parammat[,p] <- paramvec
  }
  temp <- data.frame(itnumber, populid, populWI, populBI, id, parammat)
  temp$id <- as.factor(temp$id)
  return(temp)
}


##################################################
#### calculation SD within individuals from data #
##################################################
calcSDwithin <- function(df) {
  SDwithinvec <- c()
  for (i in 1:max(as.numeric(df$id))){
    df.temp <- df[df$id==i,]
    SDwithinvec <- c(SDwithinvec, sd(df.temp[,2]))
  } 
  SDwithin <- mean(SDwithinvec)
  return(SDwithin)
}


###################################################
#### calculation SD between individuals from data #
###################################################
calcSDbetween <- function(df) {
  means <- tapply(df[,2], df[,1], mean)
  SDbetween <- sd(means)
  return(SDbetween)
}




