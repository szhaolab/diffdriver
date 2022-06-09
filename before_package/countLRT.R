library(data.table)

cmodel <- function(mutdatalist,edata){
  # i is position index, j is sample index
  lln <- function(betaf0){
    # log likelihood under null 
    delta.n <- exp(betaf0) 
    ll.n.temp <- rep(0,length(mutdatalist))
    for (t in 1:length(mutdatalist)){
      npos <- dim(mutdatalist[[t]])[1]
      ntotal <- npos * dim(mutdatalist[[t]])[2]
      nmut <- sum(mutdatalist[[t]])
      ll.n.temp[t] <- (ntotal - nmut) * log(1- BMR[t]*delta.n) + nmut * log(BMR[t]*delta.n) # bernoulli density
    }
    return(sum(ll.n.temp))
  }
  
  lla <- function(beta){
    # log likelihood under alt
    betaf0 <- beta[1]
    betag <- beta[2]
    delta.a <- exp(betaf0 + betag*edata) 
    ll.a.temp <- rep(0,length(mutdatalist))
    for (t in 1:length(mutdatalist)){
      npos <- dim(mutdatalist[[t]])[1]
      ntotal <- npos * dim(mutdatalist[[t]])[2]
      nmut_j <- colSums(mutdatalist[[t]])
      ll.a.temp[t] <- sum((npos -nmut_j) * log(1- BMR[t]*delta.a) + nmut_j * log(BMR[t]*delta.a))
    }
    return(sum(ll.a.temp))
  }
  
  resn <- optim(0, lln, method="BFGS", control=list(fnscale=-1))
  resa <- optim(c(0,0), lla, method="BFGS", control=list(fnscale=-1))
  
  teststat <- -2*(resn$value-resa$value)
  pvalue <- pchisq(teststat,df=1,lower.tail=FALSE)
  res <- list("pvalue"=pvalue, "null.beta" = resn$par, "alt.beta"=resa$par, "null.ll"= resn$value, "alt.ll" = resa$value)
  return(res)
}

cmodel.fe1 <- function(mutdatalist, edata, fehat, annodata){
  ## mutdatalist, edata, annodata should match
  ## i is position index, j is sample index
  ## use pre-fixed functional effect sizes (treat it as known)
  ## for null, also use fehat
  lln <- function(betaf0){
    ## log likelihood under null
    delta.n1 <- exp(betaf0)
    delta.n2 <- exp(betaf0 + fehat)
    # countlist <- list()
    ll.n.temp <- rep(0,length(mutdatalist))
    for (t in 1:length(mutdatalist)){
      npos1 <- dim(annodata[[t]][functypecode==7])[1]
      ntotal1 <- npos1 * dim(mutdatalist[[t]])[2]
      nmut1 <- sum(mutdatalist[[t]][1:npos1,])
      ll.n.temp1 <- (ntotal1 - nmut1) * log(1- BMR[t]*delta.n1) + nmut1 * log(BMR[t]*delta.n1) # bernoulli density
      
      npos2 <- dim(annodata[[t]][functypecode==8])[1]
      ntotal2 <- npos2 * dim(mutdatalist[[t]])[2]
      nmut2 <- sum(mutdatalist[[t]][(npos1+1):(npos1+npos2),])
      ll.n.temp2 <- (ntotal2 - nmut2) * log(1- BMR[t]*delta.n2) + nmut2 * log(BMR[t]*delta.n2) # bernoulli density
      ll.n.temp[t] <-  ll.n.temp1 + ll.n.temp2
      # countlist[[t]] <- c(npos1, npos2, nmut1, nmut2)
    }
    # print(colSums(do.call(rbind,countlist)))
    return(sum(ll.n.temp))
  }
  
  lla <- function(beta){
    ## log likelihood under alt
    betaf0 <- beta[1]
    betag <- beta[2]
    delta.a1 <- exp(betaf0 + betag*edata)
    delta.a2 <- exp(betaf0 + betag*edata + fehat)
    ll.a.temp <- rep(0, length(mutdatalist))
    # countlist <- list()
    for (t in 1:length(mutdatalist)){
      npos1 <- dim(annodata[[t]][functypecode==7])[1]
      ntotal1 <- npos1 * dim(mutdatalist[[t]])[2]
      nmut_j1 <- colSums(mutdatalist[[t]][1:npos1,])
      ll.a.temp1 <- sum((npos1 -nmut_j1) * log(1- BMR[t]*delta.a1) + nmut_j1 * log(BMR[t]*delta.a1))
      
      npos2 <- dim(annodata[[t]][functypecode==8])[1]
      ntotal2 <- npos2 * dim(mutdatalist[[t]])[2]
      nmut_j2 <- colSums(mutdatalist[[t]][(npos1+1):(npos1+npos2),])
      ll.a.temp2 <- sum((npos2 -nmut_j2) * log(1- BMR[t]*delta.a2) + nmut_j2 * log(BMR[t]*delta.a2))
      ll.a.temp[t] <- ll.a.temp1 + ll.a.temp2
      #countlist[[t]] <- c(sum(nmut_j1),sum(nmut_j2))
    }
    # print(colSums(do.call(rbind,countlist)))
    return(sum(ll.a.temp))
  }
  
  resn <- optim(0, lln, method="BFGS", control=list(fnscale=-1))
  resa <- optim(c(0,0), lla, method="BFGS", control=list(fnscale=-1))
  
  teststat<--2*(resn$value-resa$value)
  pvalue <- pchisq(teststat,df=1,lower.tail=FALSE)
  res <- list("pvalue"=pvalue, "null.beta" = resn$par, "alt.beta"=resa$par, "null.ll"= resn$value, "alt.ll" = resa$value)
  return(res)
}

cmodel.fe2 <- function(mutdatalist, edata, fehat, annodata){
  ## mutdatalist, edata, annodata should match
  ## i is position index, j is sample index
  ## use pre-fixed functional effect sizes (treat it as known)
  ## for null, use 0 for fehat
  lln <- function(betaf0){
    # log likelihood under null 
    delta.n1 <- exp(betaf0)
    delta.n2 <- exp(betaf0)
    ll.n.temp <- rep(0,length(mutdatalist))
    for (t in 1:length(mutdatalist)){
      npos1 <- dim(annodata[[t]][functypecode==7])[1]
      ntotal1 <- npos1 * dim(mutdatalist[[t]])[2]
      nmut1 <- sum(mutdatalist[[t]][1:npos1,])
      ll.n.temp1 <- (ntotal1 - nmut1) * log(1- BMR[t]*delta.n1) + nmut1 * log(BMR[t] * delta.n1) # bernoulli density
      
      npos2 <- dim(annodata[[t]][functypecode==8])[1]
      ntotal2 <- npos2 * dim(mutdatalist[[t]])[2]
      nmut2 <- sum(mutdatalist[[t]][(npos1+1):(npos1+npos2),])
      ll.n.temp2 <- (ntotal2 - nmut2) * log(1- BMR[t]*delta.n2) + nmut2 * log(BMR[t]*delta.n2) # bernoulli density
      ll.n.temp[t] <-  ll.n.temp1 + ll.n.temp2
    }
    return(sum(ll.n.temp))
  }
  
  lla <- function(beta){
    # log likelihood under alt
    betaf0 <- beta[1]
    betag <- beta[2]
    delta.a1 <- exp(betaf0 + betag*edata)
    delta.a2 <- exp(betaf0 + betag*edata +fehat) 
    ll.a.temp <- rep(0,length(mutdatalist))
    for (t in 1:length(mutdatalist)){
      npos1 <- dim(annodata[[t]][functypecode==7])[1]
      ntotal1 <- npos1 * dim(mutdatalist[[t]])[2]
      nmut_j1 <- colSums(mutdatalist[[t]][1:npos1,])
      ll.a.temp1 <- sum((npos1 -nmut_j1) * log(1- BMR[t]*delta.a1) + nmut_j1 * log(BMR[t]*delta.a1))
      
      npos2 <- dim(annodata[[t]][functypecode==8])[1]
      ntotal2 <- npos2 * dim(mutdatalist[[t]])[2]
      nmut_j2 <- colSums(mutdatalist[[t]][(npos1+1):(npos1+npos2),])
      ll.a.temp2 <- sum((npos2 -nmut_j2) * log(1- BMR[t]*delta.a2) + nmut_j2 * log(BMR[t]*delta.a2))
      ll.a.temp[t] <- ll.a.temp1 + ll.a.temp2
    }
    return(sum(ll.a.temp))
  }
  
  resn <- optim(0, lln, method="BFGS", control=list(fnscale=-1))
  resa <- optim(c(0,0), lla, method="BFGS", control=list(fnscale=-1))
  
  teststat<--2*(resn$value-resa$value)
  pvalue <- pchisq(teststat,df=1,lower.tail=FALSE)
  res <- list("pvalue"=pvalue, "null.beta" = resn$par, "alt.beta"=resa$par, "null.ll"= resn$value, "alt.ll" = resa$value)
  return(res)
}


cmodel.frac <- function(mutdatalist, edata, annodata, betaf){
  ## mutdatalist, edata, annodata should match
  ## i is position index, j is sample index
  ## allow for selection in both null and alternative, but with different fractions
  ## for selected, the effect sizes are fixed (from driverMAPS), for neutral, the effect sizes are zeros.
  ## infer using E=0 data, then test if eta is equal to the one inferred from E=0 (only works for two group comparison)
  
  betaf0 <- betaf[1]
  betaf1 <- betaf[2]
  lln <- function(eta0){
    # log likelihood under null 
    delta.n1 <- exp(eta0) * (exp(betaf0)-1) + 1 
    delta.n2 <- exp(eta0) * (exp(betaf0 + betaf1)-1) + 1
    ll.n.temp <- rep(0,length(mutdatalist))
    for (t in 1:length(mutdatalist)){
      npos1 <- dim(annodata[[t]][functypecode==7])[1]
      if (npos1 == 0) {
        ll.n.temp1 <- 0
      } else {
        ntotal1 <- npos1 * dim(mutdatalist[[t]])[2]
        nmut1 <- sum(mutdatalist[[t]][1:npos1,])
        ll.n.temp1 <- (ntotal1 - nmut1) * log(1- BMR[t]*delta.n1) + nmut1 * log(BMR[t]*delta.n1) # bernoulli density
      }

      npos2 <- dim(annodata[[t]][functypecode==8])[1]
      if (npos2 == 0) {
        ll.n.temp2 <- 0
      } else {
        ntotal2 <- npos2 * dim(mutdatalist[[t]])[2]
        nmut2 <- sum(mutdatalist[[t]][(npos1+1):(npos1+npos2),])
        ll.n.temp2 <- (ntotal2 - nmut2) * log(1- BMR[t]*delta.n2) + nmut2 * log(BMR[t]*delta.n2) # bernoulli density
      }
      ll.n.temp[t] <-  ll.n.temp1 + ll.n.temp2
    }
    return(sum(ll.n.temp))
  }
  
  lla <- function(eta){
    # log likelihood under alt
    eta0 <- eta[1]
    eta1 <- eta[2]
    delta.a1 <- (exp(eta1) * (exp(betaf0)-1) + 1) * edata + (exp(eta0) * (exp(betaf0)-1) + 1) * (1- edata)
    delta.a2 <- (exp(eta1) * (exp(betaf0 + betaf1)-1) + 1) * edata + (exp(eta0) * (exp(betaf0 + betaf1)-1) + 1) * (1- edata)
    ll.a.temp <- rep(0,length(mutdatalist))
    for (t in 1:length(mutdatalist)){
      npos1 <- dim(annodata[[t]][functypecode==7])[1]
      if (npos1 == 0) {
        ll.a.temp1 <- 0
      } else {
        ntotal1 <- npos1 * dim(mutdatalist[[t]])[2]
        nmut_j1 <- colSums(mutdatalist[[t]][1:npos1,])
        ll.a.temp1 <- sum((npos1 -nmut_j1) %*% log(1- BMR[t]*delta.a1) + nmut_j1 %*% log(BMR[t]*delta.a1))
      }
      
      npos2 <- dim(annodata[[t]][functypecode==8])[1]
      if (npos2 == 0) {
        ll.a.temp2 <- 0
      } else {
        ntotal2 <- npos2 * dim(mutdatalist[[t]])[2]
        nmut_j2 <- colSums(mutdatalist[[t]][(npos1+1):(npos1+npos2),])
        ll.a.temp2 <- sum((npos2 -nmut_j2) %*% log(1- BMR[t]*delta.a2) + nmut_j2 %*% log(BMR[t]*delta.a2))    
      }
      ll.a.temp[t] <- ll.a.temp1 + ll.a.temp2
    }
    return(sum(ll.a.temp))
  }
  resa <- optim(c(0,0), lla, method="BFGS", control=list(fnscale=-1))
  resn <- optim(0, lln, method="BFGS", control=list(fnscale=-1))

  teststat<- -2*(resn$value-resa$value)
  pvalue <- pchisq(teststat,df=1,lower.tail=FALSE)
  res <- list("pvalue"=pvalue, "null.eta0" = resn$par, "alt.eta"=resa$par, "null.ll"= resn$value, "alt.ll" = resa$value)
  return(res)
}


