#' @title gene level multiple linear regression
#' @export
mlr <- function(mut, e,covariates=1){
  # multilinear regression on gene level
  if(covariates==1){
    covariates=rep(1,length(e))
    }
  dstatus <- colSums(mut)
  dstatus[which(dstatus>0)] <- 1
  mlrfit <- lm(e ~ dstatus+covariates)
  mlrres <- summary(mlrfit)
  res <- list("res" = mlrres, "pvalue" = mlrres$coefficients[2,4])
  return(res)
}

#' @title gene level multiple linear regression, correcting for total number of mutations
#' @export
mlr.v2 <- function(mut, e, nmut,covariates=1){
  # multilinear regression adjusting for mutation rate
  # nmut is number of mutations per sample
  if(covariates==1){
    covariates=rep(1,length(e))
  }
  dstatus <- colSums(mut)
  dstatus[which(dstatus>0)] <- 1
  #dstatus <- as.factor(dstatus)
  mlrfit <- lm(e ~ dstatus + nmut +covariates)
  mlrres <- summary(mlrfit)
  res <- list("res" = mlrres, "pvalue" = mlrres$coefficients[2,4])
  return(res)
}

#' @title gene level binomial test
#' @export
genebinom <- function(mut, e){
  # binom on gene level
  mu.c <- length(e[e==1])
  mu.n <- length(e[e==0])
  gmutc <- sum(mut[,e==1])
  gmutn <- sum(mut[,e==0])
  res <- list("pvalue" = binom.test(gmutc, gmutn + gmutc, p=mu.c/(mu.c+mu.n))$p.value)
  return(res)
}

#' @title gene level logistic regression
#' @export
genelr <- function(mut, e,covariates=1){
  # logistic regression on gene level
  if(covariates==1){
    covariates=rep(1,length(e))
  }
  dstatus <- colSums(mut)
  dstatus[which(dstatus>0)] <- 1
  glrfit <- glm(e ~ dstatus+covariates, family = binomial(link = "logit"))
  glrres <- summary(glrfit)
  res <- list("res" = glrres, "pvalue" = glrres$coefficients[2,4])
  return(res)
}

#' @title gene level fisher's exact test
#' @export
genefisher <- function(mut, e){
  # fisher's test on gene level
  mu.c <- length(e[e==1])
  mu.n <- length(e[e==0])
  gmutc <- sum(mut[,e==1])
  gmutn <- sum(mut[,e==0])
  res <- list("pvalue" = fisher.test(matrix(c(gmutc, gmutn,mu.c, mu.n), ncol=2))$p.value, "count" = c(gmutc, gmutn, mu.c, mu.n))
  return(res)
}
