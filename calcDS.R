# calculation of LDA discrimination score
# imput is a dataframe with individual identity in first column and individual identity features in subsequent columns
# does not modify dataframe in any way (obtaining principal components, scaling, centering - all must be done separately)
# sets equal priors for each individual
# uses leave-one-out cross-validation

calcDS <- function (temp){
  temp[,1] <- as.factor(temp[,1])
  Y <- temp[-1]
  prior <- rep(1/nlevels(temp[,1]), nlevels(temp[,1]))
  lda.sim <- lda(Y, grouping=temp[,1], CV=T, prior=prior)
  ct <- table(temp[,1], lda.sim$class)
  result <- sum(diag(prop.table(ct)))
  return(result)
}
