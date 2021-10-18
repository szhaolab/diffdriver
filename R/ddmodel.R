#' @title diffDriver model
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
ddmodel <- function(mut, e, bmr, fe){
  
  fe <- exp(fe) #un-logs the log scale
  bmr <- exp(bmr)
  mut.t <- t(mut) #transpose of mutation matrix
  mut <- as.matrix(mut) #converts to matrix if needed?
  mut.t <-  as.matrix(mut.t)
  
  
  #Calculations completed for whole gene (using element math?)
  
  lln <- function(eta){
    # log likelihood under null
    # Under H0 (null hypothesis), i.e. there is no differential selection, α1=0. The parameter we want to estimate is α0.
    b.pi <- exp(eta)/(1+exp(eta))
    b.pi[b.pi <= 0] <- 1e-8
    b.pi[b.pi >= 1] <- 1 - 1e-8
    llmtx = log(exp(-bmr*fe)*fe*bmr*b.pi*e + exp(-bmr)*(1-b.pi*e))
    ll <- sum(llmtx)
    return(ll)
  }
  
  lla <- function(eta){
    # log likelihood under alt
    # Under H1 (alternative hypothesis), α1≠0, there are two parameters we want to estimate, α0 and α1
    eta0 <- eta[1]
    eta1 <- eta[2]
    b.pi <-  exp(eta0+eta1*e)/(1+exp(eta0+eta1*e))
    b.pi[b.pi <= 0] <- 1e-8
    b.pi[b.pi >= 1] <- 1 - 1e-8
    llmtx = log((exp(-bmr*fe)*fe*bmr*b.pi*e) + exp(-bmr)*(1-b.pi*e))
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

  resn <- optim(0, lln, method="Nelder-Mead", control=list(fnscale=-1)) # BFGS has Error: non-finite finite-difference value [2]
  resa <- optim(c(0,0), lla, method="Nelder-Mead", control=list(fnscale=-1))

  teststat<- -2*(resn$value-resa$value)
  pvalue <- pchisq(teststat,df=1,lower.tail=FALSE)
  res <- list("pvalue"=pvalue, "null.eta0" = resn$par, "alt.eta"=resa$par, "null.ll"= resn$value, "alt.ll" = resa$value)
  return(res)
}
