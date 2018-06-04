
dataz <- read.delim('./empirical datasets/CCformants.csv', sep=';', header=T)
temp <- dataz
dataz[,1] <- as.factor(dataz[,1])
dataz <- calcPCAscaled(dataz) #scale variables (zero mean and unit variance) and calculate PCA

# plot first two Principal components
plot(dataz[,2], dataz[,3], col=dataz[,1], xlab='PC1', ylab='PC2', cex=1, cex.lab=1)

# calculate univariate identity metrics
calcF(dataz)
calcPICbetweenmeans(dataz) # gives nonsense values; needs to be calculated on non-transformed data
calcPICbetweenmeans(temp)
calcHSnpergroup(dataz,sum=F)

# calculate multivariate identity metrics
calcDS(dataz)
calcHSnpergroup(dataz, sum=T)
calcMI(dataz)
calcHM(dataz)


