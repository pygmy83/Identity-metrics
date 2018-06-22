calcMI <- function(temp){
  temp[,1] <- as.factor(temp[,1])
  Y <- temp[-1]
  prior <- rep(1/nlevels(temp[,1]), nlevels(temp[,1]))
  lda.sim <- lda(Y, grouping=temp[,1], CV=T, prior=prior)
  MI <- mutinformation(temp[,1], lda.sim$class) #in nats
  MI <- natstobits(MI) #conversion to bits
  return(MI)
}
