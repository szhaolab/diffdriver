


#' get pi
#' @param alpha
#' @param e
#' @return
#' @export
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

#' Title
#'
#' @param b
#' @param zpost
#' @param rate.s0
#' @param ll.n
#' @param mutidx
#'
#' @return
#' @export
#'
#' @examples
q_pos_nl <- function(b, zpost, rate.s0, ll.n, mutidx){
  ll.s <- get_ll_s_nl(b, rate.s0, mutidx)
  q <- sum(zpost[ ,1] * ll.s + zpost[ ,2] * ll.n)
  return(q)
}

#' Title
#'
#' @param p
#' @param rate.s0
#' @param ll.n
#' @param mutidx
#'
#' @return
#' @export
#'
#' @examples
dd_loglik_nl <- function(p, rate.s0, ll.n, mutidx,e){
  beta0 <- p[1]
  alpha <- p[2:3]
  pi <- get_pi_nl(alpha, e)
  ll.s <- get_ll_s_nl(beta0, rate.s0, mutidx)
  ll <- sum(log(pi * exp(ll.s) + (1-pi) * exp(ll.n)))
}

#' Title
#'
#' @param p
#' @param rate.n
#' @param rate.s0
#' @param ll.n
#' @param mutidx
#' @param type
#'
#' @return
#' @export
#'
#' @examples
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
#  reslg <- nnet::nnet.default(lg.x, zpost, size = 0,
#                        skip = TRUE, softmax = TRUE, censored = FALSE,
#                        rang = 0, trace=FALSE)
#  coef <- coefficients(reslg) #  b->o1    i1->o1     b->o2    i1->o2
#  # 1: intercept for category 1, 2: slope for variable 1 in category 1.
#  # 3: intercept for category 2, 2: slope for variable 1 in category 2.
#  alpha <- c(coef[1] - coef[3],coef[2] - coef[4])
#  if (type == "null"){
#    alpha <- c(sum(alpha), 0)
#  }
#
#  pnew <- c(beta0, alpha)
#

 alpha <- brglm(zpost~ lg.x,family="binomial")$coefficients
 if (type == "null"){
     alpha <- c(alpha[1], 0)
   }
   pnew <- c(beta0, alpha)
   return(pnew)(pnew)
}

#' Title
#'
#' @param beta0
#' @param alpha
#' @param rate.n
#' @param rate.s0
#' @param ll.n
#' @param mutidx
#' @param type
#' @param maxit
#' @param tol
#'
#' @return
#' @export
#'
#' @examples
#dd_EM_ordinary <- function(beta0 = 0, alpha = c(0,0), rate.n, rate.s0, ll.n, mutidx, type = c("null", "alt"), maxit = 100, tol = 1e-3,mut,e){
#  ll_rec <- rep(0, maxit)
#  p_rec <- NULL
#
#  # initialize
#  p <- c(beta0, alpha)
#
#  for (i in 1:maxit){
#    ll <- dd_loglik(p, rate.s0, ll.n, mutidx,e)
#    ll_rec[i] <- ll
#    cat("iteration ", i,"; loglikelihood:", ll, "\n")
#    pnew <- dd_EM_update(p, rate.n, rate.s0, ll.n, mutidx, type = type,mut=mut,e=e)
#    p_rec <- rbind(p_rec, pnew)
#
#    if  (dist(rbind(pnew, p)) < tol){
#      break
#    }
#    p <- pnew
#  }
#
#  return(list("loglikelihood" = ll, "beta0" = p[1], "alpha" = p[2:3], "loglik_rec" = ll_rec, "param_rec" = p_rec))
#}
#

#' Title
#'
#' @param beta0
#' @param alpha
#' @param rate.n
#' @param rate.s0
#' @param ll.n
#' @param mutidx
#' @param type
#' @param maxit
#' @param tol
#'
#' @return
#' @export
#'
#' @examples
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
#' @export
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

