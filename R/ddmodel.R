#' @title diffDriver model
#' @description This model is applied on data of a single gene. It will infer effect size for both sample-level variable and positional level functional annotations. We used an EM algorithm to infer parameters.
#' @param mut a matrix of mutation status 0 or 1, rows positions, columns are samples.
#' @param e a vector,phenotype of each sample,
#'  should match the columns of \code{mut} and \code{mr}
#' @param mr a matrix, mutation rate of each sample at each mutation (log scale) that is not dependent on sample level factor
#' @param fe, a vector, increased mutation rate at each position, depending on e (log scale),
#'  should match the rows of \code{mut} and \code{mr}
#' @export
ddmodel <- function(mut, e, mr, fe){
  rate.n <- as.matrix(exp(mr))
  rate.s0 <- as.matrix(exp(fe) * rate.n)
  ll.n <- colSums(log(rate.n * mut +  (1-rate.n) * (1-mut)))

  get_ll_s <- function(b){
    rate.s <- rate.s0 * exp(b)
    # log likelihood for each sample under selection
    colSums(log(rate.s * mut +  (1-rate.s) * (1-mut)))
  }

  get_pi_null <- function(alpha, e){
    pi <- exp(alpha)/(1 + exp(alpha))
    rep(pi, length(e))
  }

  get_pi_alt <- function(alpha, e){
    alpha0 <- alpha[1]
    alpha1 <- alpha[2]
    exp(alpha0 + alpha1 * e)/(1 + exp(alpha0 + alpha1 * e))
  }

  q_pos <- function(b, zpost){
    ll.s <- get_ll_s(b)
    q <- sum(zpost * ll.s + (1-zpost) * ll.n)
    return(q)
  }

  ddmodel_EM <- function(niter = 30, type = c("null", "alt")){
    # Initialize
    beta0 <- 0

    if (type == "null"){
      alpha = 0
      pi <- get_pi_null(alpha, e)
    }

    if (type == "alt"){
      alpha = c(0,0)
      pi <- get_pi_alt(alpha, e)
    }

    for (i in 1:niter){
      # update z_i
      ll.s <- get_ll_s(beta0)
      zpost <- (pi * exp(ll.s) + (1-pi) * exp(ll.n))/sum(pi * exp(ll.s) + (1-pi) * exp(ll.n))

      # update beta0
      res <- optim(0, q_pos, zpost = zpost, method = "BFGS", control=list(fnscale=-1))
      beta0 <- res$par

      # update alpha
      if (type == "null"){
        reslg <- multinom(pi ~ 1)
      }
      if (type == "alt"){
        reslg <- multinom(pi ~ e)
      }
      alpha <- coefficients(reslg)
      pi <- reslg$fitted.values
    }

    ll <- sum(pi * exp(ll.s) + (1-pi) * exp(ll.n))
    return(list("pi" = pi, "loglikelihood" = ll, "beta0" = beta0, "alpha" = alpha))
  }

  res.null <- ddmodel_EM(type = "null")
  res.alt <- ddmodel_EM(type = "alt")
  teststat<- -2*(res.null$loglikelihood - res.alt$loglikelihood)
  pvalue <- pchisq(teststat, df=1, lower.tail=FALSE)
  res <- list("pvalue"=pvalue, "res.null" = res.null, "res.alt"=res.alt)
  return(res)

}


#' @title diffDriver model with effect size for positional functional annotations fixed
#' @description This function uses the model as cmodel.frac,
#'  but generalizes to take more than 1 functional categories.
#'  This model is applied on data of a single gene
#' @param mut a matrix of mutation status 0 or 1, rows positions, columns are samples.
#' @param e a vector,phenotype of each sample,
#'  should match the columns of \code{mut} and \code{mr}
#' @param mr a matrix, mutation rate of each sample at each mutation (log scale) that is not dependent on sample level factor
#' @param fe, a vector, increased mutation rate at each position, depending on e (log scale),
#'  should match the rows of \code{mut} and \code{mr}
#' @export
ddmodel_fixb <- function(mut, e, mr, fe){
  mr <- exp(mr)
  fe <- exp(fe)

  rate.s <- as.matrix(fe * mr)
  rate.n <-  as.matrix(mr)
  rate.s[rate.s <= 0] <- 0
  rate.s[rate.s >= 1] <- 1
  rate.n[rate.n <= 0] <- 0
  rate.n[rate.n >= 1] <- 1
  ll.s <- log(rate.s * mut +  (1-rate.s) * (1-mut)) # matrix, log likelihood for each (i,j) under selection
  ll.n <- log(rate.n * mut +  (1-rate.n) * (1-mut)) # matrix, log likelihood for each (i,j) no selection
  l.s <- exp(colSums(ll.s))
  l.n <- exp(colSums(ll.n))

  lln <- function(param){
    # log likelihood under null
    # Under H0 (null hypothesis), i.e. there is no differential selection, α1=0. The parameter we want to estimate is α0.
    alpha <- param[1]
    pi1 <- exp(alpha)/(1 + exp(alpha))
    ll = sum(log(l.s * pi1 + l.n * (1-pi1)))
    return(ll)
  }

  lla <- function(param){
    # log likelihood under alt
    # Under H1 (alternative hypothesis), α1≠0, there are two parameters we want to estimate, α0 and α1

    alpha0 <- param[1]
    alpha1 <- param[2]
    pi1 <-  exp(alpha0 + alpha1 * e)/(1 + exp(alpha0 + alpha1 * e))
    ll = sum(log(l.s * pi1 + l.n * (1-pi1)))
    return(ll)
  }

  resn <- optim(c(0), lln, method="Nelder-Mead", control=list(fnscale=-1)) # BFGS has Error: non-finite finite-difference value [2]
  resa <- optim(c(0,0), lla, method="Nelder-Mead", control=list(fnscale=-1))

  teststat<- -2*(resn$value-resa$value)
  pvalue <- pchisq(teststat,df=1,lower.tail=FALSE)
  res <- list("pvalue"=pvalue, "null.param" = resn$par, "alt.param"=resa$par, "null.ll"= resn$value, "alt.ll" = resa$value)
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
    ll <- sum(llmtx)
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
    ll <- sum(llmtx)
    return(ll)
  }

  resn <- optim(0, lln, method="Nelder-Mead", control=list(fnscale=-1)) # BFGS has Error: non-finite finite-difference value [2]
  resa <- optim(c(0,0), lla, method="Nelder-Mead", control=list(fnscale=-1))

  teststat<- -2*(resn$value-resa$value)
  pvalue <- pchisq(teststat,df=1,lower.tail=FALSE)
  res <- list("pvalue"=pvalue, "null.eta0" = resn$par, "alt.eta"=resa$par, "null.ll"= resn$value, "alt.ll" = resa$value)
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

  lln <- function(eta){
    # log likelihood under null
    b.pi <- exp(eta) * (fe - 1) * bmr[,1] + bmr[,1]
    b.pi[b.pi <= 0] <- 1e-8
    b.pi[b.pi >= 1] <- 1 - 1e-8

    llmtx <- npos * log(1 - b.pi) + mutpos * log(b.pi)
    ll <- sum(llmtx)
    return(ll)
  }

  lla <- function(eta){
    # log likelihood under alt
    eta0 <- eta[1]
    eta1 <- eta[2]
    b.pi0 <- t(exp(eta0) * (fe - 1) * bmr[,1] + bmr[,1])
    b.pi1 <- t(exp(eta1) * (fe - 1) * bmr[,1] + bmr[,1])
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
