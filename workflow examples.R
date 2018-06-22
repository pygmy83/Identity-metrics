
dataz <- read.delim('./empirical datasets/CCformants.csv', sep=';', header=T) # load dataset
temp <- dataz # backup of original values
dataz[,1] <- as.factor(dataz[,1]) # first column needs to be factor with individual identity labels
dataz <- calcPCAscaled(dataz) #scale variables (zero mean and unit variance) and calculate PCA
plot(dataz[,2], dataz[,3], col=dataz[,1], xlab='PC1', ylab='PC2', cex=1, cex.lab=1) # illustration of variation in dataset based on first two Principal components

# calculate univariate identity metrics
calcF(dataz)
calcPICbetweenmeans(dataz) # gives nonsense values; PIC cannot handle values around zero properly; needs to be calculated on non-transformed data
calcPICbetweenmeans(temp) # using original, non-transformed data for PIC
calcHSnpergroup(dataz,sum=F) # gives HS for each variable separately

# calculate multivariate identity metrics
calcDS(dataz)
calcHSnpergroup(dataz, sum=T) # gives total HS for all variables or for only those that significantly differ between individuals (in one-way ANOVA with identity as independent factor)
calcMI(dataz)
calcHM(dataz)

# calculate summary for multiple datasets
destination <- './empirical datasets examples' # the folder with multiple CSV files each representing a single dataset  
DatasetSummary(destination, pca=T)
#write.table(temp, paste0(destination, '/datasetsummary.csv'), sep=';')


# predict HS and DS values based on loess models
nindividuals <- 40 # number of individuals in the analysis - current models require values betwee 5 and 40
nobservations <- 20 # number of observations per individual in the analysis - current models require values betwee 5 and 20
DS <- 0.80 # total score from discrimination analysis - must be between 0 and 1
datazpred <- data.frame(w=nindividuals, x=nobservations, y=DS)
predict(DStoHSloess, datazpred)

nindividuals <- 40 # number of individuals in the analysis - current models require values betwee 5 and 40
nobservations <- 20 # number of observations per individual in the analysis - current models require values betwee 5 and 20
HS <- 12.9 # HS value; current model require values between 0 and 12.9
datazpred <- data.frame(w=nindividuals, x=nobservations, z=HS)
predict(HStoDSloess, datazpred) 



###########################################################################################
# Calculate univariate metrics on simulated datasets for a given combination of parameters#
###########################################################################################

varmean <- 1000 # mean value for a variable
idRange <- c(1) # range of values for individuality
iRange <- c(5,10,30) # range of values for number of individuals
oRange <- c(10) # range of values for number of observations
itRange <- c(1:10) # range of values for number of iterations
result <- ItterateDatasetsUni(idRange, iRange, oRange, it, varmean)

y <- result$HSngroups
x <- result$i
plotmeans(y~x)

# settings for the analysis documenting relationships between metrics and parameters
# varmean <- 1000 # mean value for a variable
# idRange <- c(0.01, 1, 2.5, 5, 10) # range of values for individuality
# iRange <- c(5, 10, 15, 20, 25, 30, 35, 40) # range of values for number of individuals
# oRange <- c(4, 8, 12, 16, 20) # range of values for number of observations per individual
# itRange <- c(1:20) # range of values for number of iterations


##############################################################################################
# Calculate multivariate metrics on simulated datasets for a given combinations of parameters#
##############################################################################################

iRange <- c(5,10,30) # range of values for number of individuals
oRange <- c(5) # range of values for number of observations per individual 
pRange <- c(2) # range of values for number of variables describing individual identity = pcs
covRange <- c(0) # range of values for covariance between variables
idRange <- c(1) # range of values for individuality in data
itRange <- c(1:10)# range of values for number of iterations
result <- ItterateDatasets(covRange,idRange,iRange,oRange,pRange,itRange)

y <- result$DS
x <- result$i
plotmeans(y~x)

# settings for the analysis documenting relationships between metrics and parameters
# iRange <- c(5, 10, 15, 20, 25, 30, 35, 40) 
# oRange <- c(4, 8, 12, 16, 20) #, 12, 16, 20) 
# pRange <- c(2, 4, 6, 8, 10) #, 6, 8, 10)  
# itRange <- c(1:20)
# idRange <- c(0.01, 1, 2.5, 5, 10) 
# covRange <- c(0, 0.25, 0.5, 0.75, 1) 

# analysis to get HS values up to about 10 for DS and HS conversions
# iRange <- c(5, 15, 20, 25, 30, 35, 40) # range for number of individuals:
# oRange <- c(4, 8, 12, 16, 20) # range for number of observations = calls
# pRange <- c(2, 4, 6, 8, 10) #, 6, 8, 10)  # range for number of parameters = pcs
# itRange <- c(1:10)# range for number of iterations
# idRange <- c(0.1, 0.25, 0.5, 0.75, 1, 1.33, 1.66, 2)
# covRange <- c(0)

