calcMI <- function(temp){
  temp[,1] <- as.factor(temp[,1])
  Y <- temp[-1]
  prior <- rep(1/nlevels(temp[,1]), nlevels(temp[,1]))
  lda.sim <- lda(Y, grouping=temp[,1], CV=T, prior=prior)
  MI <- mutinformation(temp[,1], lda.sim$class) #in nats
  MI <- natstobits(MI) #conversion to bits
  return(MI)
}

calcMI2 <- function(temp){
  temp[,1] <- as.factor(temp[,1])
  Y <- temp[-1]
  prior <- rep(1/nlevels(temp[,1]), nlevels(temp[,1]))
  lda.sim <- lda(Y, grouping=temp[,1], CV=T, prior=prior)
  ct <- table(temp[,1], lda.sim$class)
  MI <- f_rez(temp[,1], lda.sim$class)
  return(MI) #first number is NMI second number is AMI
}



# functions used to calculate MI
# http://stackoverflow.com/questions/21831953/r-package-available-for-adjusted-mutual-information

f_rez <- function(v1,v2){
  s1 <- tabulate(v1);
  s2 <- tabulate(v2);
  l1 <- length(s1)
  l2 <- length(s2)
  N <- length(v1)
  tij <- f_nij(v1,v2,l1,l2)   #contingency table n(i,j)=t(i,j). this would be equivalent with table(v1,v2)
  mi <- mutinformation(v1,v2) #function for Mutual Information from package infotheo
  h1 <- -sum(s1*log(s1/N))/N
  h2 <- -sum(s2*log((s2+1)/(N+l1)))/N  # CHANGED!!! originally: h2 <- -sum(s2*log(s2/N))/N problem that s2 is sometimes 0########################!!!!
  nmi <- mi/max(h1,h2)        # NMI Normalized MI
  ami <- NA
  #emi <- f_emi(s1,s2,l1,l2,N) # EMI Expected MI
  #ami <- (mi-emi)/max(h1,h2)  #AMI Adjusted MI
  return(c(nmi,ami,mi)) ### changed to return MI as well ##################################################!!!!  
}

f_nij <- function(v1,v2,l1,l2){     #contingency table n(i,j)=t(i,j)
  m <- matrix(0,l1,l2)
  for (i in 1:length(v1)){
    m[v1[i],v2[i]] <- m[v1[i],v2[i]] +1
  }
  m
}

f_emi <- function(s1,s2,l1,l2,n){    #expected mutual information
  s_emi <- 0
  for(i in 1:l1){
    for (j in 1:l2){
      min_nij <- max(1,s1[i]+s2[j]-n)
      max_nij <- min(s1[i],s2[j])
      n.ij <- seq(min_nij, max_nij)   #sequence of consecutive numbers
      t1<- (n.ij / n) * log((n.ij * n) / (s1[i]*s2[j]))
      t2 <- exp(lfactorial(s1[i]) + lfactorial(s2[j]) + lfactorial(n - s1[i]) + lfactorial(n - s2[j]) - lfactorial(n) - lfactorial(n.ij) - lfactorial(s1[i] - n.ij) - lfactorial(s2[j] - n.ij) - lfactorial(n - s1[i] - s2[j] + n.ij))
      emi <- sum(t1*t2)
      s_emi <- s_emi + emi
    }
  }
  return(s_emi)
}
