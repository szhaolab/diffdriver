#' @title diffDriver model
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
ddmodel <- function(mut, e, mr, fe){
  mr <- exp(mr)
  fe <- exp(fe) #un-logs the log scale

  #Calculations completed for whole gene (using element math?)

  lln <- function(param){
    # log likelihood under null
    # Under H0 (null hypothesis), i.e. there is no differential selection, α1=0. The parameter we want to estimate is α0.
    alpha <- param[1]
    #fe0 <- param[2]

    rate.s <-  fe * mr
    rate.n <-  mr

    rate.s[rate.s <= 0] <- 1e-8
    rate.s[rate.s >= 1] <- 1 - 1e-8
    rate.n[rate.n <= 0] <- 1e-8
    rate.n[rate.n >= 1] <- 1 - 1e-8

    l.s <- rate.s * mut +  (1-rate.s) * (1-mut) #matrix, likelihood for each (i,j) under selection
    l.n <- rate.n * mut +  (1-rate.n) * (1-mut) #matrix, likelihood for each (i,j) no selection

    l.s <- t(l.s)
    l.n <- t(l.n)

    pi1 <- exp(alpha)/(1 + exp(alpha))
    llmtx = log(l.s * pi1 + l.n * (1-pi1))
    ll <- sum(llmtx)
    return(ll)
  }

  lla <- function(param){
    # log likelihood under alt
    # Under H1 (alternative hypothesis), α1≠0, there are two parameters we want to estimate, α0 and α1

    alpha0 <- param[1]
    alpha1 <- param[2]
    # fe0 <- param[3]

    rate.s <-  fe * mr
    rate.n <-  mr

    rate.s[rate.s <= 0] <- 1e-8
    rate.s[rate.s >= 1] <- 1 - 1e-8
    rate.n[rate.n <= 0] <- 1e-8
    rate.n[rate.n >= 1] <- 1 - 1e-8

    l.s <- rate.s * mut +  (1-rate.s) * (1-mut) #matrix, likelihood for each (i,j) under selection
    l.n <- rate.n * mut +  (1-rate.n) * (1-mut) #matrix, likelihood for each (i,j) no selection

    l.s <- t(l.s)
    l.n <- t(l.n)

    pi1 <-  exp(alpha0 + alpha1 * e)/(1+exp(alpha0 + alpha1 * e))
    llmtx = log(l.s * pi1 + l.n * (1-pi1))
    ll <- sum(llmtx)
    return(ll)
  }

  resn <- optim(c(0), lln, method="Nelder-Mead", control=list(fnscale=-1)) # BFGS has Error: non-finite finite-difference value [2]
  resa <- optim(c(0,0), lla, method="Nelder-Mead", control=list(fnscale=-1))

  teststat<- -2*(resn$value-resa$value)
  pvalue <- pchisq(teststat,df=1,lower.tail=FALSE)
  res <- list("pvalue"=pvalue, "null.param" = resn$par, "alt.param"=resa$par, "null.ll"= resn$value, "alt.ll" = resa$value)
  return(res)
}
