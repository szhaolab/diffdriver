
get_pi <- function(alpha, e){
  alpha0 <- alpha[1]
  alpha1 <- alpha[2]
  if (alpha0>10 | alpha1>10){
    return(p=1)
  }
  else{
    p=exp(alpha0 + alpha1 * e)/(1 + exp(alpha0 + alpha1 * e))
  }
}


q_pos <- function(b, zpost, rate.s0, mut,mutidx,ll.n){
  ll.s <- get_ll_s(b,mut,mutidx, rate.s0)
  q <- sum(zpost[ ,1] * ll.s + zpost[ ,2] * ll.n)
  return(q)
}

dd_loglik <- function(p, rate.s0, ll.n,mutidx, mut, e){
  beta0 <- p[1]
  alpha <- p[2:3]
  pi <- get_pi(alpha, e)
  ll.s <- get_ll_s(beta0, mut,mutidx, rate.s0)
  ll <- sum(log(pi * exp(ll.s) + (1-pi) * exp(ll.n)))
}


get_ll_s <- function(b, mut, mutidx, rate_s0){
  rate_s <- rate_s0 * exp(b)
  rmtx <- -rate_s
  rmtx[mutidx] <- mut[mutidx]*log(rate_s[mutidx]) - rate_s[mutidx]
  # log likelihood for each sample under selection
  colSums(rmtx)
}


dd_EM_update <- function(p, rate.n, rate.s0, ll.n, mutidx,type = c("null", "alt"),mut,e){
  # p: beta0, alpha
  beta0 <- p[1]
  alpha <- p[2:3]

  # update z_i
  ll.s <- get_ll_s(b=beta0,mut=mut,mutidx=mutidx, rate_s0=rate.s0)
  pi <- get_pi(alpha=alpha, e=e)
  zpost <- cbind(pi * exp(ll.s) , (1-pi) * exp(ll.n)) # 1st column selection, 2nd column neutral
  zpost <- zpost/rowSums(zpost)

  # update beta0
  #bi=(nrow(mutidx) - sum(rate.n %*% zpost[,2,drop=F]))/sum(rate.s0 %*% zpost[,1,drop=F])
  #beta0.init <- ifelse( bi>0, log(bi),rnorm(1))
  beta0.init <- 1
  res <- optim(beta0.init, q_pos, zpost = zpost, rate.s0 = rate.s0, mut=mut,ll.n = ll.n, mutidx=mutidx, control=list(fnscale=-1))
  beta0 <- res$par

  # update alpha
  if (type == "null"){
    lg.x <- rep(1, ncol(mut))
  }
  if (type == "alt"){
    lg.x <- e
  }

  alpha <- brglm::brglm(zpost~ lg.x,family="binomial")$coefficients

  if (type == "null"){
    alpha <- c(alpha[1], 0)
  }

  pnew <- c(beta0, alpha)

  return(pnew)
}


dd_EM_ordinary <- function(beta0 = 0, alpha = c(0,0), rate.n, rate.s0, ll.n, mutidx, type = c("null", "alt"), maxit = 100, tol = 1e-3,mut,e){
  ll_rec <- rep(0, maxit)
  p_rec <- NULL

  # initialize
  p <- c(beta0, alpha)

  for (i in 1:maxit){
    ll <- dd_loglik(p, rate.s0, ll.n, mutidx,e)
    ll_rec[i] <- ll
    cat("iteration ", i,"; loglikelihood:", ll, "\n")
    pnew <- dd_EM_update(p, rate.n, rate.s0, ll.n, mutidx, type = type,mut=mut,e=e)
    p_rec <- rbind(p_rec, pnew)

    if  (dist(rbind(pnew, p)) < tol){
      break
    }
    p <- pnew
  }

  return(list("loglikelihood" = ll, "beta0" = p[1], "alpha" = p[2:3], "loglik_rec" = ll_rec, "param_rec" = p_rec))
}



dd_squarEM <- function(beta0 = 0, alpha = c(0,0), rate.n, rate.s0, ll.n,mutidx, type = c("null", "alt"), mut,e, maxit = 100, tol = 1e-3){

  # initialize
 p <- c(beta0, alpha)

  # EM
  res <- SQUAREM::squarem(p=p, rate.n = rate.n, rate.s0=rate.s0, ll.n=ll.n, mutidx=mutidx, type = type,mut=mut,e=e, fixptfn=dd_EM_update, control=list(tol=tol, maxiter = maxit))
  p <- res$par
  ll <- dd_loglik(p=p, rate.s0=rate.s0, ll.n=ll.n, mutidx=mutidx,mut=mut,e=e)

  return(list("loglikelihood" = ll, "beta0" = p[1], "alpha" = p[2:3]))
}

#' @title diffDriver model
#' @description This model is applied on data of a single gene. It will infer effect size for both sample-level variable and positional level functional annotations. We used an EM algorithm to infer parameters.
#' @param mut a matrix of mutation status 0 or 1, rows positions, columns are samples.
#' @param e a vector,phenotype of each sample,
#'  should match the columns of \code{mut} and \code{mr}
#' @param mr a matrix, mutation rate of each sample at each mutation (log scale) that is not dependent on sample level factor
#' @param fe, a vector, increased mutation rate at each position, depending on e (log scale),
#'  should match the rows of \code{mut} and \code{mr}
#' @export
ddmodel <- function(mut, e, mr, fe, label, ...){
  rate.n <- exp(mr)
  rate.s0 <- exp(fe) * rate.n

  ## aggregate the duplicate rows
  rate.n <- aggregate(rate.n,by=list(label),sum,na.rm=T)[,-c(1)]
  rate.s0 <- aggregate(rate.s0,by=list(label),sum,na.rm=T)[,-c(1)]
  mut <- aggregate(as.matrix(mut),by=list(label),sum,na.rm=T)[,-c(1)]
  mutidx <- which(mut!=0, arr.ind = T)

  ## Parameter estimation based on Poisson likelihood

  ll.n <- colSums(mut*log(rate.n)-rate.n)

  res.null <- dd_squarEM(rate.n = rate.n, rate.s0 = rate.s0, ll.n=ll.n,mutidx=mutidx, type = "null",mut=mut,e=e, ...)

  res.alt <- dd_squarEM(rate.n = rate.n, rate.s0 = rate.s0, ll.n=ll.n,mutidx=mutidx, type = "alt",mut=mut,e=e, ...)

  ## Get p value

  teststat<- -2*(res.null$loglikelihood - res.alt$loglikelihood)
  pvalue <- pchisq(teststat, df=1, lower.tail=FALSE)
  res <- list("pvalue"=pvalue, "res.null" = res.null, "res.alt"=res.alt)
  return(res)
}



#' @title diffDriver model only for binary phenotype.
#' @description This function uses the model as cmodel.frac,
#'  but generalizes to take more than 1 functional categories.
#'  This model is applied on data of a single gene
#' @param mut a matrix of mutation status 0 or 1
#' @param e a vector,phenotype of each sample,
#'  should match the columns of \code{mut} and \code{bmr}
#' @param bmr a matrix, background mutation rate of each sample at each mutation (log scale)
#' @param fe, a vector, increased mutation rate at each mutation, due to functional effect (log scale),
#'  should match the rows of \code{mut} and \code{bmr}
#' @export
ddmodel_binary <- function(mut, e, bmr, fe){

  fe <- exp(fe)
  bmr <- exp(bmr)
  mut.t <- t(mut)
  mut <- as.matrix(mut)
  mut.t <-  as.matrix(mut.t)

  lln <- function(eta){
    # log likelihood under null
    b.pi <- exp(eta) * (fe - 1) * bmr + bmr
    b.pi[b.pi <= 0] <- 1e-8
    b.pi[b.pi >= 1] <- 1 - 1e-8
    llmtx <- (1- mut) * log(1 - b.pi) + mut * log(b.pi)
    ll <- sum(llmtx,na.rm=T)
    return(ll)
  }

  lla <- function(eta){
    # log likelihood under alt
    eta0 <- eta[1]
    eta1 <- eta[2]
    b.pi0 <- t(exp(eta0) * (fe - 1) * bmr + bmr)
    b.pi1 <- t(exp(eta1) * (fe - 1) * bmr + bmr)
    b.pi0[b.pi0 <= 0] <- 1e-8
    b.pi1[b.pi1 <= 0] <- 1e-8
    b.pi0[b.pi0 >= 1] <- 1 - 1e-8
    b.pi1[b.pi1 >= 1] <- 1 - 1e-8
    b.pi <-  b.pi0 + e * (b.pi1- b.pi0)
    llmtx <- (1- mut.t) * log(1- b.pi) + mut.t * log(b.pi)
    ll <- sum(llmtx,na.rm=T)
    return(ll)
  }
  resn <- optim(0, lln, method="Nelder-Mead", control=list(fnscale=-1)) # BFGS has Error: non-finite finite-difference value [2]
  resa <- optim(c(0,0), lla, method="Nelder-Mead", control=list(fnscale=-1))
  probnull <- mean(exp(resn$par)*(fe-1)*bmr+bmr)
  probalt1 <- mean(exp(resa$par[1])*(fe-1)*bmr+bmr)
  probalt2 <- mean(exp(resa$par[2])*(fe-1)*bmr+bmr)

  teststat<- -2*(resn$value-resa$value)
  pvalue <- pchisq(teststat,df=1,lower.tail=FALSE)
all<- c(pvalue,resn$par,resa$par,resn$value,resa$value)
names(all)=c("pvalue","null_para","e=0","e=1","null_likelihood","alt_likelihood")
  res <- list("pvalue"=pvalue, "null.eta0" = resn$par, "alt.eta"=resa$par, "null.ll"= resn$value, "alt.ll" = resa$value,"all"=all)
  return(res)
}


#' @title diffDriver model only for binary phenotype, assuming bmr the same across samples.
#' @description This function uses the model as cmodel.frac,
#'  but generalizes to take more than 1 functional categories.
#'  This model is applied on data of a single gene. This should give the same results as ddmodel_binary defined above.
#' @param mut a matrix of mutation status 0 or 1
#' @param e a vector,phenotype of each sample,
#'  should match the columns of \code{mut} and \code{bmr}
#' @param bmr a matrix, background mutation rate of each sample at each mutation (log scale), as we assume bmr the same across samples, only the first column will be used.
#' @param fe, a vector, increased mutation rate at each mutation, due to functional effect (log scale),
#'  should match the rows of \code{mut} and \code{bmr}
#' @export
ddmodel_binary_simple <- function(mut, e, bmr, fe){

  fe <- exp(fe)
  bmr <- exp(bmr)
  mut <- as.matrix(mut)

  mutpos <- rowSums(mut)
  npos <- ncol(mut) - mutpos
  mutpos.e0 <- rowSums(mut[, e==0])
  npos.e0 <- ncol(mut[, e==0]) -  mutpos.e0
  mutpos.e1 <- rowSums(mut[, e==1])
  npos.e1 <- ncol(mut[, e==1]) - mutpos.e1
index0=min(which(e==0))
index1=min(which(e==1))


  lln <- function(eta){
    # log likelihood under null
    b.pi0 <- exp(eta) * (fe - 1) * bmr[,index0] + bmr[,index0]
    b.pi1 <- exp(eta) * (fe - 1) * bmr[,index1] + bmr[,index1]


    b.pi0[b.pi0 <= 0] <- 1e-8
    b.pi0[b.pi0 >= 1] <- 1 - 1e-8
    b.pi1[b.pi1 <= 0] <- 1e-8
    b.pi1[b.pi1 >= 1] <- 1 - 1e-8

    llmtx0 <- npos.e0 * log(1- b.pi0) + mutpos.e0 * log(b.pi0)
    llmtx1 <- npos.e1 * log(1- b.pi1) + mutpos.e1 * log(b.pi1)
    ll <- sum(llmtx0,na.rm=T)+sum(llmtx1,na.rm=T)
    return(ll)
  }


  lla <- function(eta){
    # log likelihood under alt
    eta0 <- eta[1]
    eta1 <- eta[2]
    index0=min(which(e==0))
    index1=min(which(e==1))
    b.pi0 <- t(exp(eta0) * (fe - 1) * bmr[,index0] + bmr[,index0])
    b.pi1 <- t(exp(eta1) * (fe - 1) * bmr[,index1] + bmr[,index1])
    b.pi0[b.pi0 <= 0] <- 1e-8
    b.pi1[b.pi1 <= 0] <- 1e-8
    b.pi0[b.pi0 >= 1] <- 1 - 1e-8
    b.pi1[b.pi1 >= 1] <- 1 - 1e-8
    llmtx0 <- npos.e0 * log(1- b.pi0) + mutpos.e0 * log(b.pi0)
    llmtx1 <- npos.e1 * log(1- b.pi1) + mutpos.e1 * log(b.pi1)
    ll <- sum(llmtx0) + sum(llmtx1)
    return(ll)
  }

  resn <- optim(0, lln, method="Nelder-Mead", control=list(fnscale=-1)) # BFGS has Error: non-finite finite-difference value [2]
  resa <- optim(c(0,0), lla, method="Nelder-Mead", control=list(fnscale=-1))

  teststat<- -2*(resn$value-resa$value)
  pvalue <- pchisq(teststat,df=1,lower.tail=FALSE)
  res <- list("pvalue"=pvalue, "null.eta0" = resn$par, "alt.eta"=resa$par, "null.ll"= resn$value, "alt.ll" = resa$value)
  return(res)
}


# note this version does not estimate beta. It doesn't use EM, but use optim to estimate alpha.
ddmodel_conti_simple <- function(mut, e, bmr, fe){

  fe <- exp(fe)
  bmr <- exp(bmr)
  mut <- as.matrix(mut)
  mut1=data.frame()
  mut0=data.frame()
  aa=data.frame()
  ubmr=unique(bmr)[,1]
  n=ncol(bmr)

  for (i in l:length(ubmr)){
    index=which(bmr[,1]==ubmr[i])
    ife=fe[index]
    ufe=unique(ife)
    for (j in 1:length(ufe)){
      index=which(bmr[,1]==ubmr[i] & fe==ufe[j] )
      imut=mut[index,]
      mut1=rbind(mut1,rowSums(imut))
      mut0=rbind(mut0,nrow(imut)-rowSums(imut))
      aa=rbind(aa,c(ubmr[i],ubmr[i]*ufe[j]))
    }
  }

  lln <- function(eta){
    # log likelihood under null
    b.pi0 <- aa[,1]^mut1*(1-aa[,1])^mut0
    b.pi1 <- aa[,2]^mut1*(1-aa[,2])^mut0


    b.pi0[b.pi0 <= 0] <- 1e-8
    b.pi0[b.pi0 >= 1] <- 1 - 1e-8
    b.pi1[b.pi1 <= 0] <- 1e-8
    b.pi1[b.pi1 >= 1] <- 1 - 1e-8
    pi= exp(eta)/(1+exp(eta))
    ll <- sum(log(b.pi0+b.pi1),na.rm=T)
    return(ll)
  }


  lla <- function(eta){
    # log likelihood under alt
    eta0 <- eta[1]
    eta1 <- eta[2]
    index0=min(which(e==0))
    index1=min(which(e==1))
    b.pi0 <- t(exp(eta0) * (fe - 1) * bmr[,index0] + bmr[,index0])
    b.pi1 <- t(exp(eta1) * (fe - 1) * bmr[,index1] + bmr[,index1])
    b.pi0[b.pi0 <= 0] <- 1e-8
    b.pi1[b.pi1 <= 0] <- 1e-8
    b.pi0[b.pi0 >= 1] <- 1 - 1e-8
    b.pi1[b.pi1 >= 1] <- 1 - 1e-8
    llmtx0 <- npos.e0 * log(1- b.pi0) + mutpos.e0 * log(b.pi0)
    llmtx1 <- npos.e1 * log(1- b.pi1) + mutpos.e1 * log(b.pi1)
    ll <- sum(llmtx0,na.rm=T) + sum(llmtx1,na.rm=T)
    return(ll)
  }

  resn <- optim(0, lln, method="Nelder-Mead", control=list(fnscale=-1)) # BFGS has Error: non-finite finite-difference value [2]
  resa <- optim(c(0,0), lla, method="Nelder-Mead", control=list(fnscale=-1))

  teststat<- -2*(resn$value-resa$value)
  pvalue <- pchisq(teststat,df=1,lower.tail=FALSE)
  res <- list("pvalue"=pvalue, "null.eta0" = resn$par, "alt.eta"=resa$par, "null.ll"= resn$value, "alt.ll" = resa$value)
  return(res)
}


