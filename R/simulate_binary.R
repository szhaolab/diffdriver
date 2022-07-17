#' @title simulate phenotype and mutations for one gene, using one functional annotation.
#' @description simulate mutations for one gene. The sample level phenotype are binary (two groups.).
#' Using only one functional category for position level annotation (missense and loss of function).
#' @param sgdata a list, each list item for one nttype
#' @param bmrpars a vector, each item is bmr in log scale for one nttype
#' @param betaf0, shift of mutation rate from BMR, log scale, shared in all samples.
#' @param Nsample, total number of samples
#' @param beta_gc, effect size for functional covariates, log scale. Right now we only have one functional covariate, that is whether the mutation is missense or loss of function mutation. beta_gc[1] indicates the shift of mutation rate for missense ( coded as 7 in functype code column in `sgdata`), beta_gc[2] indicates the shift of mutation rate for loss of function ( coded as 8 in functype code column in `sgdata`). log scale.
#' @param para par[1] is the fraction of positively selected samples in group E=1; par[2] is the
#' fraction of negatively selected samples in group E=0
#' @param hotspot hotspot[1] is the probability of being hotspot for a given position. hotspot[2] is the log size
#' for hotspots.    
#' @import Matrix data.table
#' @export
simulate_1funcv <- function(sgdata, bmrpars, betaf0, Nsample, beta_gc,para,hotspot){
  Nsamplec <- round(Nsample/2) # number of samples with phenotype E=1 (the rest will be 0)
  Nsamplen <- Nsample-Nsamplec
  edata <- c(rep(1,Nsamplec),rep(0,Nsamplen))
  
  Nsamplec.ps <- rbinom(1, Nsamplec, para[1]) # number of hotspots position in positively selected samples  in group E=1
  #Nsamplec.ps.hot <- rbinom(1, Nsamplec.ps, hotspot[1]) # number of hotspots position in positively selected samples  in group E=1
  #Nsamplec.ps.reg=Nsamplec.ps-Nsample.ps.hot
  Nsamplen.ps <- rbinom(1, Nsamplen, para[2]) # number of positively selected samples in group E=0
  Nsample.ps <- Nsamplec.ps + Nsamplen.ps
  Nsample.neu <- Nsample - Nsample.ps
  mutlist <- list()
  countlist <- list()
  annodata <- list()
  bmrmtxlist <- list()
  hotsize <- c()
  for (t in 1:length(sgdata)){ # Simulate mutation data. t: nucleotide change type
    tnpos1 <- dim(sgdata[[t]][functypecode==7])[1]
    tnpos2 <- dim(sgdata[[t]][functypecode==8])[1]
    k1=round(tnpos1*hotspot[1])
    k2=tnpos1-k1
    k3=round(tnpos2*hotspot[1])
    k4=tnpos2-k3
    annodata[[t]] <- rbind(sgdata[[t]][functypecode==7], sgdata[[t]][functypecode==8])
    
    
    mutc1.hot <- rsparsematrix(k1, Nsample.ps, nnz=rbinom(1, Nsample.ps * k1, exp(bmrpars[t])*exp(betaf0)*exp(beta_gc[1]))*exp(hotspot[2]), rand.x=NULL)
    mutc1.reg <- rsparsematrix(k2, Nsample.ps, nnz=rbinom(1, Nsample.ps * k2, exp(bmrpars[t])*exp(betaf0)*exp(beta_gc[1])), rand.x=NULL)
    mutc1=rbind(mutc1.hot,mutc1.reg)
    
    mutc2.hot <- rsparsematrix(k3, Nsample.ps, nnz=rbinom(1, Nsample.ps * k3, exp(bmrpars[t])*exp(betaf0)*exp(beta_gc[1] + beta_gc[2])*exp(hotspot[2])), rand.x=NULL)
    mutc2.reg <- rsparsematrix(k4, Nsample.ps, nnz=rbinom(1, Nsample.ps * k4, exp(bmrpars[t])*exp(betaf0)*exp(beta_gc[1] + beta_gc[2])), rand.x=NULL)
    mutc2=rbind(mutc2.hot,mutc2.reg)
    
    mutc <- rbind(mutc1, mutc2)

    
    mutn1.hot <- rsparsematrix(k1, Nsample.neu, nnz=rbinom(1, Nsample.neu * k1, exp(bmrpars[t])*exp(betaf0))*exp(hotspot[2]), rand.x=NULL)
    mutn1.reg <- rsparsematrix(k2, Nsample.neu, nnz=rbinom(1, Nsample.neu * k2, exp(bmrpars[t])*exp(betaf0)), rand.x=NULL)
    mutn1=rbind(mutn1.hot,mutn1.reg)
  
    
    mutn2.hot <- rsparsematrix(k3, Nsample.neu, nnz=rbinom(1, Nsample.neu * k3, exp(bmrpars[t])*exp(betaf0)*exp(hotspot[2])), rand.x=NULL)
    mutn2.reg <- rsparsematrix(k4, Nsample.neu, nnz=rbinom(1, Nsample.neu * k4, exp(bmrpars[t])*exp(betaf0)), rand.x=NULL)
    mutn2=rbind(mutn2.hot,mutn2.reg)
    
    mutn <- rbind(mutn1, mutn2)

    mutc.out <- cbind(mutc[, 1:Nsamplec.ps], mutn[, 1:(Nsamplec-Nsamplec.ps)])
    mutn.out <- cbind(mutc[, (Nsamplec.ps + 1):Nsample.ps], mutn[, (Nsamplec-Nsamplec.ps+1):Nsample.neu])

    mutlist[[t]] <- cbind(mutc.out,mutn.out)
    countlist[[t]] <- c(tnpos1, tnpos2,sum(mutc1),sum(mutc2),sum(mutn1), sum(mutn2), sum(mutc.out), sum(mutn.out), sum(mutc.out[1:tnpos1,]), sum(mutn.out[1:tnpos1,]))
    bmrmtxlist[[t]] <- matrix(bmrpars[t], ncol = ncol(mutlist[[t]]), nrow = nrow(mutlist[[t]]))
    aa=hotspot[2]*c(rep(1,floor(tnpos1*hotspot[1])),rep(0,floor(tnpos1*(1-hotspot[1]))),
     rep(1,floor(tnpos2*hotspot[1])),rep(0,floor(tnpos1*(1-hotspot[1])))) 
    hotsize=c(hotsize,aa)
      }

  avbetaf1 <- log(exp(beta_gc[1]) * Nsample.ps/Nsample + Nsample.neu/Nsample)
  avbetaf2 <- log(exp(beta_gc[1] + beta_gc[2]) * Nsample.ps/Nsample + Nsample.neu/Nsample)
  pos1pos2ratio <- colSums(do.call(rbind, countlist))[1]/colSums(do.call(rbind, countlist))[2]
  avbetaf1f2 <- log((pos1pos2ratio*exp(avbetaf1) + exp(avbetaf2))/(pos1pos2ratio+1))
  betaf1f2 <- log((pos1pos2ratio*exp(beta_gc[1]) + exp(beta_gc[2]))/(pos1pos2ratio+1))
  simdata <- list("mutlist"= mutlist, "pheno" = edata, "annodata" = annodata, "bmrpars" = bmrpars, "bmrmtxlist" = bmrmtxlist, "para" = para, "hotsize"=hotsize, "efsize" = list( "betaf0" = betaf0,  "beta_gc" = beta_gc, "avbetaf1" = avbetaf1, "avbetaf2" = avbetaf2, "avbetaf1f2" = avbetaf1f2,"betaf1f2"=betaf1f2))
  return(simdata)
}
