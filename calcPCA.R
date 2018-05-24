calcPCA <- function (temp){
  id <- temp[,1]
  pcscores <- prcomp(temp[,-1], center=F, scale=F, tol=0)$x
  temp <- data.frame(id,pcscores)
  return(temp)
}


calcPCAscaled <- function (temp){
  id <- temp[,1]
  pcscores <- prcomp(temp[,-1], center=T, scale=T, tol=0)$x
  temp <- data.frame(id,pcscores)
  return(temp)
}
