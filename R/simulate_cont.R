#' @title simulate phenotype and mutations for one gene, using one functional annotation.
#' @description simulate mutations for one gene. The sample level phenotype are binary (two groups.).
#' Using only one functional category for position level annotation (missense and loss of function).
#' @param sgdata a list, each list item for one nttype
#' @param bmrpars a vector, each item is bmr in log scale for one nttype
#' @param betaf0, shift of mutation rate from BMR, log scale, shared in all samples.
#' @param Nsample, total number of samples
#' @param beta_gc, effect size for functional covariates, log scale. Right now we only have one functional covariate, that is whether the mutation is missense or loss of function mutation. beta_gc[1] indicates the shift of mutation rate for missense ( coded as 7 in functype code column in `sgdata`), beta_gc[2] indicates the shift of mutation rate for loss of function ( coded as 8 in functype code column in `sgdata`). log scale.
#' @param par, the parameters for generating the data,  in which par[1] is the mean of phenotype, par[2] the sd of phenotype,
#' par[3] the intercept of logistic regression and par[4] the slope of the logistic regression. Here the outcome of logistic regression is the
#' positive sample and the independent variable is the phenotype.
#' @param hotspot hotspot[1] is the probability of being hotspot for a given position. hotspot[2] is the log size
#' for hotspots.
#' @import Matrix data.table
#' @export
simulate_2funcv <- function(binary=F,sgdata, bmrpars, betaf0=2, Nsample, beta_gc, para,hotseq, hmm){
  if (binary==T){
    Nsamplec <- round(Nsample/2) # number of samples with phenotype E=1 (the rest will be 0)
    Nsamplen <- Nsample-Nsamplec
    phenotype <- c(rep(1,Nsamplec),rep(0,Nsamplen))
    ss=ifelse(phenotype==1,sample(c(0,1),size=Nsamplec,replace=T,prob = c(1-para[1],para[1])),
              sample(c(0,1),size=Nsamplen,replace=T,prob = c(1-para[2],para[2])))

  }else{
  phenotype=rnorm(Nsample,mean = para[1],sd=para[2])
  pp=exp(para[3]+para[4]*phenotype)/(1+exp(para[3]+para[4]*phenotype))
  ss=ifelse(runif(Nsample)-pp<0,1,0)
  }

  index=which(ss==1)
  phenotype=c(phenotype[index],phenotype[-index])
  Nsample.ps=sum(ss)
  Nsample.neu <- Nsample - Nsample.ps





  mutlist <- list()
  countlist <- list()
  annodata <- list()
  bmrmtxlist <- list()
  hotsize <- list()
  for (t in 1:length(sgdata)) {


    ssgdata=merge(sgdata[[t]],hotseq,by="start")

    tnpos1 <- dim(sgdata[[t]][functypecode==7])[1]
    tnpos2 <- dim(sgdata[[t]][functypecode==8])[1]

    seqt=ssgdata$seqt

    k1=sum(seqt[ssgdata$functypecode==7])
    k2=tnpos1-k1
    k3=sum(seqt[ssgdata$functypecode==8])
    k4=tnpos2-k3

    pp2=exp(bmrpars[t])*exp(betaf0)*exp(beta_gc[1])
    pp3=exp(bmrpars[t])*exp(betaf0)*exp(beta_gc[1] + beta_gc[2])
    pp1=exp(bmrpars[t])*exp(betaf0)
    hotpp1= min(pp1*exp(hmm[9]),1)
    hotpp2= min(pp2*exp(hmm[9]),1)
    hotpp3= min(pp3*exp(hmm[9]),1)

    annodata[[t]] <- rbind(sgdata[[t]][functypecode==7], sgdata[[t]][functypecode==8])
    mut1=matrix(nrow = tnpos1,ncol = Nsample.ps)
    mut2=matrix(nrow = tnpos2,ncol = Nsample.ps)
    mut3=matrix(nrow = tnpos1,ncol = Nsample.neu)
    mut4=matrix(nrow = tnpos2,ncol = Nsample.neu)

    for (j in 1:Nsample.ps) {
      if (k1>0){
        mut1[1:k1,j]= sample(c(0,1),k1, replace=T,prob= c(1-hotpp2,hotpp2))
        mut1[(1+k1):tnpos1,j]= sample(c(0,1),k2, replace = T, prob = c(1-pp2,pp2))
      }else{
        mut1[,j]=sample(c(0,1),k2, replace = T, prob = c(1-pp2,pp2))
      }

      if (k3>0){
        mut2[1:k3,j]= sample(c(0,1),k3, replace=T,prob= c(1-hotpp3,hotpp3))
        mut2[(1+k3):tnpos2,j]= sample(c(0,1),k4, replace = T, prob = c(1-pp3,pp3))
      }else{
      mut2[,j]=sample(c(0,1),k4, replace = T, prob = c(1-pp3,pp3))
      }
    }
    mut.ps=rbind(mut1,mut2)


    for (j in 1:Nsample.neu) {
      mut3[,j]= sample(c(0,1),tnpos1, replace=T,prob= c(1-pp1,pp1))
      mut4[,j]= sample(c(0,1),tnpos2, replace=T,prob= c(1-pp1,pp1))
    }
   mut.neu=rbind(mut3,mut4)

    mutlist[[t]] <- cbind(mut.ps,mut.neu)
    countlist[[t]] <- c(tnpos1, tnpos2,sum(mut.ps),sum(mut.neu))
    bmrmtxlist[[t]] <- matrix(bmrpars[t], ncol = ncol(mutlist[[t]]), nrow = nrow(mutlist[[t]]))
    hotsize[[t]]=hmm[9]*c(rep(1,k1),rep(0,k2),
                    rep(1,k3),rep(0,k4))
    }

  avbetaf1 <- log(exp(beta_gc[1]) * Nsample.ps/Nsample + Nsample.neu/Nsample)
  avbetaf2 <- log(exp(beta_gc[1] + beta_gc[2]) * Nsample.ps/Nsample + Nsample.neu/Nsample)
  pos1pos2ratio <- colSums(do.call(rbind, countlist))[1]/colSums(do.call(rbind, countlist))[2]
  avbetaf1f2 <- log((pos1pos2ratio*exp(avbetaf1) + exp(avbetaf2))/(pos1pos2ratio+1))
  betaf1f2 <- log((pos1pos2ratio*exp(beta_gc[1]) + exp(beta_gc[2]))/(pos1pos2ratio+1))
  simdata <- list("mutlist"= mutlist, "hotspots"=hotsize, "pheno" = phenotype, "annodata" = annodata, "bmrpars" = bmrpars, "bmrmtxlist" = bmrmtxlist, "para"=para, "hotsize"=hotsize,"phenotype"=phenotype,  "efsize" = list( "betaf0" = betaf0,  "beta_gc" = beta_gc, "avbetaf1" = avbetaf1, "avbetaf2" = avbetaf2, "avbetaf1f2" = avbetaf1f2,"betaf1f2"=betaf1f2),"nsample"=c(Nsample.ps,Nsample.neu))
  return(simdata)
}



#' revised simulation function
#'
#' @param binary
#' @param sgdata
#' @param bmrpars
#' @param betaf0
#' @param Nsample
#' @param beta_gc
#' @param para
#' @param hotseq
#' @param hmm
#'
#' @return
#' @export
#'
#' @examples
simulate_1funcv <- function(binary=F,sgdata, bmrpars, betaf0=2, Nsample, beta_gc, para,hot=0, hmm){
  if (binary==T){
    Nsamplec <- round(Nsample/2) # number of samples with phenotype E=1 (the rest will be 0)
    Nsamplen <- Nsample-Nsamplec
    phenotype <- c(rep(1,Nsamplec),rep(0,Nsamplen))
    ss=ifelse(phenotype==1,sample(c(0,1),size=Nsamplec,replace=T,prob = c(1-para[1],para[1])),
              sample(c(0,1),size=Nsamplen,replace=T,prob = c(1-para[2],para[2])))

  }else{
    phenotype=rnorm(Nsample,mean = para[1],sd=para[2])
    pp=exp(para[3]+para[4]*phenotype)/(1+exp(para[3]+para[4]*phenotype))
    ss=ifelse(runif(Nsample)-pp<0,1,0)
  }

  index=which(ss==1)
  phenotype=c(phenotype[index],phenotype[-index])
  Nsample.ps=sum(ss)
  Nsample.neu <- Nsample - Nsample.ps


hotseq=hotspotseq(hmm,sgdata)
if (hot==0){
hotseq[,2]=hotseq[,2]*0
}

mutpslist <- list()
mutneulist <- list()
  mutlist <- list()
  countlist <- list()
  annodata <- list()
  bmrmtxlist <- list()
  hotsize <- list()
  for (t in 1:length(sgdata)) {
    ssgdata=merge(sgdata[[t]],hotseq,by="start")

    tnpos1 <- dim(sgdata[[t]][functypecode==7])[1]
    tnpos2 <- dim(sgdata[[t]][functypecode==8])[1]

    seqt=ssgdata$seqt

    k1=sum(seqt[ssgdata$functypecode==7])
    k2=tnpos1-k1
    k3=sum(seqt[ssgdata$functypecode==8])
    k4=tnpos2-k3

    pp2=exp(bmrpars[t])*exp(betaf0)*exp(beta_gc[1])
    pp3=exp(bmrpars[t])*exp(betaf0)*exp(beta_gc[1] + beta_gc[2])
    pp1=exp(bmrpars[t])*exp(betaf0)
    hotpp= min(pp1*exp(hmm[9]),1)

    annodata[[t]] <- rbind(sgdata[[t]][functypecode==7], sgdata[[t]][functypecode==8])


      if (k1>0){
        mut11= rsparsematrix(k1,Nsample.ps,nnz = rbinom(1, Nsample.ps * k1, hotpp),rand.x=NULL)
        mut12= rsparsematrix(k2,Nsample.ps,nnz = rbinom(1, Nsample.ps * k2, pp2),rand.x=NULL)
        mut1=rbind(mut11,mut12)
      }else{
        mut1= rsparsematrix(k2,Nsample.ps,nnz = rbinom(1, Nsample.ps * k2, pp2),rand.x=NULL)
      }

      if (k3>0){
        mut21= rsparsematrix(k3,Nsample.ps,nnz = rbinom(1, Nsample.ps * k3, hotpp),rand.x=NULL)
        mut22= rsparsematrix(k4,Nsample.ps,nnz = rbinom(1, Nsample.ps * k4, pp3),rand.x=NULL)
        mut2=rbind(mut21,mut22)
      }else{
        mut2=  rsparsematrix(k4,Nsample.ps,nnz = rbinom(1, Nsample.ps * k4, pp3),rand.x=NULL)
      }

    mut.ps=rbind(mut1,mut2)


 mut.neu= rsparsematrix(tnpos1+tnpos2,Nsample.neu,nnz = rbinom(1, Nsample.neu * (tnpos1+tnpos2), pp1),rand.x=NULL)

mutpslist[[t]] <- mut.ps
mutneulist[[t]] <- mut.neu
    mutlist[[t]] <- cbind(mut.ps,mut.neu)
    countlist[[t]] <- c(tnpos1, tnpos2,sum(mut.ps),sum(mut.neu))
    bmrmtxlist[[t]] <- matrix(bmrpars[t], ncol = ncol(mutlist[[t]]), nrow = nrow(mutlist[[t]]))
    hotsize[[t]]=hmm[9]*c(rep(1,k1),rep(0,k2),
                          rep(1,k3),rep(0,k4))
  }

  avbetaf1 <- log(exp(beta_gc[1]) * Nsample.ps/Nsample + Nsample.neu/Nsample)
  avbetaf2 <- log(exp(beta_gc[1] + beta_gc[2]) * Nsample.ps/Nsample + Nsample.neu/Nsample)
  pos1pos2ratio <- colSums(do.call(rbind, countlist))[1]/colSums(do.call(rbind, countlist))[2]
  avbetaf1f2 <- log((pos1pos2ratio*exp(avbetaf1) + exp(avbetaf2))/(pos1pos2ratio+1))
  betaf1f2 <- log((pos1pos2ratio*exp(beta_gc[1]) + exp(beta_gc[2]))/(pos1pos2ratio+1))
  simdata <- list("mutpslist"=mutpslist,"mutneulist"=mutneulist,"mutlist"= mutlist, "hotseq"=hotseq,"hotsize"=hotsize, "pheno" = phenotype, "annodata" = annodata, "bmrpars" = bmrpars, "bmrmtxlist" = bmrmtxlist, "para"=para, "efsize" = list( "betaf0" = betaf0,  "beta_gc" = beta_gc, "avbetaf1" = avbetaf1, "avbetaf2" = avbetaf2, "avbetaf1f2" = avbetaf1f2,"betaf1f2"=betaf1f2),"nsample"=c(Nsample.ps,Nsample.neu))
  return(simdata)
}


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
simulate_3funcv <- function(sgdata, bmrpars, betaf0=0.5, Nsample=1000, beta_gc=c(0,1), fracc=0.8, fracn=0.2){
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

  simdata <- list("mutlist"= mutlist, "pheno" = edata, "annodata" = annodata, "bmrpars" = bmrpars, "bmrmtxlist" = bmrmtxlist, "fracc" = fracc, "fracn" = fracn, "efsize" = list( "betaf0" = betaf0,  "beta_gc" = beta_gc, "avbetaf1" = avbetaf1, "avbetaf2" = avbetaf2, "avbetaf1f2" = avbetaf1f2),"nsample"=c(Nsample.ps,Nsample.neu))
  return(simdata)
}

