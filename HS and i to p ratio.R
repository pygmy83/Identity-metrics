# Setting for analysis investigating 
itopRange <- c(0.5, 1, 1.5, 2, 3, 5) # 'i to p ratio' is number of individuals (i) divided by number of variables (p)
iRange <- c(3, 5, 8, 10, 15, 20, 25, 30, 40, 50, 60, 100 )
covRange <- c(0, 0.25, 0.5, 0.75, 1)
idRange <- c(0.01, 1, 2.5, 5, 10)
pRange <- c(5,10,20)
itRange <- c(1:100)
dataz <- NULL

# calculate HS for simulated datasets
set.seed(1)
tottime <- proc.time()
for (itop in 1:length(itopRange)){
    for (it in itRange) {
      for (p in pRange) {
        i <- ceiling(p*itopRange[itop])#o <- sample (oRange, 1)
        o <- 10
        cov <- sample (covRange, 1)
        id <- sample (idRange, 1)
        temp <- GenerateMultivariateDataset(i, o, p, cov, id)
        temp <- calcPCA(temp)
        HS <- calcHSnpergroup(temp)[2]
        dataz <- rbind(dataz, data.frame(cov, id,i, o, p, it, itopRange[itop], HS))
      }
    }
  #}
}
print('total time:'); print(proc.time()-tottime)

names(dataz)[1:7] <- c('cov', 'id', 'i', 'o', 'p', 'it', 'itop')
#write.table(dataz, 'data-itop.csv', sep=';', row.names=F)


# graph and results for the analysis
dataz <- read.delim('data-itop-1000iterations.csv', header=T, sep=';')
temp <- dataz
y <- temp$HS
x <- temp$itop
plotmeans(y~x, n.label=F, bars=T, barcol='black', ylim=c(0,20), 
          ylab='HS', xlab='i to p ratio (N individuals / N parameters)',
          main='100 iterations')
abline(v=2, lty = 2, col='grey')
text(x=0.5, y=19, cex=0.8, labels = 'Fewer individuals \nthan variables', pos=4, offset = 0)
arrows(1.8, 17, x1 = 0.7, y1 = 17, length = 0.15, angle = 30, code = 2, col = 'black')


summary(aov(y~as.factor(x), data=temp))
TukeyHSD(aov(y~as.factor(x), data=temp))
