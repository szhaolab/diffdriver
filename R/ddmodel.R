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
