#### function definition ###
GenerateUnivariate <- function(nindivs,ncalls,betweenM, individuality) {
  
  id <- c()
  paramvec <- c()
  means <- rnorm(nindivs, betweenM, 1)
  indcodes <- seq(from=1, to=nindivs, by=1)
  indsds <- means # just to build same matrix as for means
  indsds <- rep(sd(means)/individuality, nindivs) # calculates sd of means for each variable and sets the individual sd to be individuality times the sd between means
  
  
  for (i in 1:nindivs) {
    paramvec <- c(paramvec, rnorm(ncalls, means[i], indsds[i]))
    id <- c(id, rep(indcodes[i], ncalls))
  }
  
  temp <- data.frame(id, paramvec) # pair the matrix with individual codes
  temp$id <- as.factor(temp$id) # individual is factor
  return(temp)
  
}


#### function test ####
nindivs <- 5 # the dataset will consist of 5 individuals
ncalls <- 10 # there will be 10 observations per individual 
betweenM <- 500 # the variable mean = 500
individuality <- 2 # between-individual SD will be twice larger than within-individual SD
temp <- GenerateUnivariate(nindivs, ncalls, betweenM, individuality)
temp
boxplot(temp$paramvec~temp$id, xlab = 'Individual no.', ylab = 'variable 1', main = paste0('Individuality = ', individuality))

rm (nindivs, ncalls, betweenM, individuality, temp)


