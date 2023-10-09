
# this version of ddmodel does not aggregate mutations with the same annotation.
# it does not take the label argument.
# it gives the same p value, but the likelihood is different.
# ddmodel is based on poisson likelihood, this one is based on bernoulli.
# ddmodel is faster (2x) than dddmodel_nl.

get_pi_nl <- function(alpha, e){
  alpha0 <- alpha[1]
  alpha1 <- alpha[2]
  if (alpha0>10 | alpha1>10){
    return(p=1)
  }
  else{
    p=exp(alpha0 + alpha1 * e)/(1 + exp(alpha0 + alpha1 * e))
  }
}

q_pos_nl <- function(b, zpost, rate.s0, ll.n, mutidx){
  ll.s <- get_ll_s_nl(b, rate.s0, mutidx)
  q <- sum(zpost[ ,1] * ll.s + zpost[ ,2] * ll.n)
  return(q)
}


get_ll_s_nl <- function(b, rate_s0, mutidx){
  rate_s <- rate_s0 * exp(b)
  rmtx <- log(1-rate_s)
  rmtx[mutidx] <- log(rate_s[mutidx])
  # log likelihood for each sample under selection
  colSums(rmtx) # faster than `colSums(log(rate.s * mut +  (1-rate.s) * (1-mut)))`
}

dd_loglik_nl <- function(p, rate.s0, ll.n, mutidx,e){
  beta0 <- p[1]
  alpha <- p[2:3]
  pi <- get_pi_nl(alpha, e)
  ll.s <- get_ll_s_nl(beta0, rate.s0, mutidx)
  ll <- sum(log(pi * exp(ll.s) + (1-pi) * exp(ll.n)))
}


dd_EM_update_nl <- function(p, rate.n, rate.s0, ll.n, mutidx, type = c("null", "alt"),mut,e){
  # p: beta0, alpha
  beta0 <- p[1]
  alpha <- p[2:3]

  # update z_i
  ll.s <- get_ll_s_nl(beta0, rate.s0, mutidx)
  pi <- get_pi_nl(alpha, e)
  zpost <- cbind(pi * exp(ll.s) , (1-pi) * exp(ll.n)) # 1st column selection, 2nd column neutral
  zpost <- zpost/rowSums(zpost)

  # update beta0
  bi=(nrow(mutidx) - sum(rate.n %*% zpost[,2,drop=F]))/sum(rate.s0 %*% zpost[,1,drop=F])
  beta0.init <- ifelse( bi>0, log(bi),rnorm(1))
  res <- optim(beta0.init, q_pos_nl, zpost = zpost, rate.s0 = rate.s0, ll.n = ll.n, mutidx = mutidx, method = "BFGS", control=list(fnscale=-1))
  beta0 <- res$par

  # update alpha
  if (type == "null"){
    lg.x <- rep(1, ncol(mut))
  }
  if (type == "alt"){
    lg.x <- e
  }

  alpha <- brglm(zpost~ lg.x,family="binomial")$coefficients
  if (type == "null"){
     alpha <- c(alpha[1], 0)
   }
  pnew <- c(beta0, alpha)
  return(pnew)
}


dd_squarEM_nl <- function(beta0 = 0, alpha = c(0,0), rate.n, rate.s0, ll.n, mutidx, type = c("null", "alt"), mut,e, maxit = 100, tol = 1e-3){

  # initialize
  p <- c(beta0, alpha)
  # EM
  res <- SQUAREM::squarem(p=p, rate.n = rate.n, rate.s0=rate.s0, ll.n=ll.n, mutidx=mutidx, type = type,mut=mut,e=e, fixptfn=dd_EM_update_nl, control=list(tol=tol, maxiter = maxit))
  p <- res$par
  ll <- dd_loglik_nl(p, rate.s0, ll.n, mutidx,e)

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
#' @noRd
ddmodel_nl <- function(mut, e, mr, fe, ...){
  rate.n <- as.matrix(exp(mr))
  rate.s0 <- as.matrix(exp(fe) * rate.n)
  ll.n <- colSums(log(rate.n * mut +  (1-rate.n) * (1-mut)))
  mutidx <- which(mut!=0, arr.ind = T)

  # res.null <- dd_EM_ordinary(rate.n = rate.n, rate.s0 = rate.s0, ll.n=ll.n, mutidx=mutidx, type = "null", ...)
  res.null <- dd_squarEM_nl(rate.n = rate.n, rate.s0 = rate.s0, ll.n=ll.n, mutidx=mutidx, type = "null",mut=mut,e=e, ...)
  # res.alt <- dd_EM_ordinary(rate.n = rate.n, rate.s0 = rate.s0, ll.n=ll.n, mutidx=mutidx, type = "alt",  ...)
  res.alt <- dd_squarEM_nl(rate.n = rate.n, rate.s0 = rate.s0, ll.n=ll.n, mutidx=mutidx, type = "alt",mut=mut,e=e, ...)
  teststat<- -2*(res.null$loglikelihood - res.alt$loglikelihood)
  pvalue <- pchisq(teststat, df=1, lower.tail=FALSE)
  res <- list("pvalue"=pvalue, "res.null" = res.null, "res.alt"=res.alt)
  return(res)
}

