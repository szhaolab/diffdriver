# same as March 2021 version

#' @title simulate phenotype and mutations for one gene, using one functional annotation.
#' @description simulate mutations for one gene. The sample level phenotype are binary (two groups.).
#' Using only one functional category for position level annotation (missense and loss of function).
#' @param sgdata a list, each list item for one nttype
#' @param bmrpars a vector, each item is bmr in log scale for one nttype
#' @param betaf0, shift of mutation rate from BMR, log scale, shared in all samples.
#' @param Nsample, total number of samples
#' @param beta_gc, effect size for functional covariates, log scale. Right now we only have one functional covariate, that is whether the mutation is missense or loss of function mutation. beta_gc[1] indicates the shift of mutation rate for missense ( coded as 7 in functype code column in `sgdata`), beta_gc[2] indicates the shift of mutation rate for loss of function ( coded as 8 in functype code column in `sgdata`). log scale.
#' @param fracc, fraction of positively selected samples in group E=1
#' @param fracn, fraction of positively selected samples in group E=0
#' @import Matrix data.table
#' @export
simulate_1funcv <- function(sgdata, bmrpars, betaf0=0.5, Nsample=1000, beta_gc=c(0,1), fracc=0.8, fracn=0.2){
  Nsamplec <- round(Nsample/2) # number of samples with phenotype E=1 (the rest will be 0)
  Nsamplen <- Nsample-Nsamplec
  edata <- c(rep(1,Nsamplec),rep(0,Nsamplen))

  Nsamplec.ps <- rbinom(1, Nsamplec, fracc) # number of positively selected samples in group E=1
  Nsamplen.ps <- rbinom(1, Nsamplen, fracn) # number of positively selected samples in group E=0
  Nsample.ps <- Nsamplec.ps + Nsamplen.ps
  Nsample.neu <- Nsample - Nsample.ps

  mutlist <- list()
  countlist <- list()
  annodata <- list()
  bmrmtxlist <- list()
  for (t in 1:length(sgdata)){ # Simulate mutation data. t: nucleotide change type
    tnpos1 <- dim(sgdata[[t]][functypecode==7])[1]
    tnpos2 <- dim(sgdata[[t]][functypecode==8])[1]
    annodata[[t]] <- rbind(sgdata[[t]][functypecode==7], sgdata[[t]][functypecode==8])
    mutc1 <- rsparsematrix(tnpos1, Nsample.ps, nnz=rbinom(1, Nsample.ps * tnpos1, exp(bmrpars[t])*exp(betaf0)*exp(beta_gc[1])), rand.x=NULL)
    mutc2 <- rsparsematrix(tnpos2, Nsample.ps, nnz=rbinom(1, Nsample.ps * tnpos2, exp(bmrpars[t])*exp(betaf0)*exp(beta_gc[1] + beta_gc[2])), rand.x=NULL)
    mutc <- rbind(mutc1, mutc2)
    mutn1 <- rsparsematrix(tnpos1, Nsample.neu, nnz=rbinom(1, Nsample.neu * tnpos1, exp(bmrpars[t])*exp(betaf0)), rand.x=NULL)
    mutn2 <- rsparsematrix(tnpos2, Nsample.neu, nnz=rbinom(1, Nsample.neu * tnpos2, exp(bmrpars[t])*exp(betaf0)), rand.x=NULL)
    mutn <- rbind(mutn1, mutn2)


    mutc.out <- cbind(mutc[, 1:Nsamplec.ps], mutn[, 1:(Nsamplec-Nsamplec.ps)])
    mutn.out <- cbind(mutc[, (Nsamplec.ps + 1):Nsample.ps], mutn[, (Nsamplec-Nsamplec.ps+1):Nsample.neu])

    mutlist[[t]] <- cbind(mutc.out,mutn.out)
    countlist[[t]] <- c(tnpos1, tnpos2,sum(mutc1),sum(mutc2),sum(mutn1), sum(mutn2), sum(mutc.out), sum(mutn.out), sum(mutc.out[1:tnpos1,]), sum(mutn.out[1:tnpos1,]))
    bmrmtxlist[[t]] <- matrix(bmrpars[t], ncol = ncol(mutlist[[t]]), nrow = nrow(mutlist[[t]]))
  }

  avbetaf1 <- log(exp(beta_gc[1]) * Nsample.ps/Nsample + Nsample.neu/Nsample)
  avbetaf2 <- log(exp(beta_gc[1] + beta_gc[2]) * Nsample.ps/Nsample + Nsample.neu/Nsample)
  pos1pos2ratio <- colSums(do.call(rbind, countlist))[1]/colSums(do.call(rbind, countlist))[2]
  avbetaf1f2 <- log((pos1pos2ratio * exp(avbetaf1) + exp(avbetaf2))/(pos1pos2ratio+1))

  simdata <- list("mutlist"= mutlist, "pheno" = edata, "annodata" = annodata, "bmrpars" = bmrpars, "bmrmtxlist" = bmrmtxlist, "fracc" = fracc, "fracn" = fracn, "efsize" = list( "betaf0" = betaf0,  "beta_gc" = beta_gc, "avbetaf1" = avbetaf1, "avbetaf2" = avbetaf2, "avbetaf1f2" = avbetaf1f2))
  return(simdata)
}
