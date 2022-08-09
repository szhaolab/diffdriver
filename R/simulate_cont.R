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
simulate_2funcv <- function(sgdata, bmrpars, betaf0=2, Nsample, beta_gc, para,hotseq, hmm){
  phenotype=rnorm(Nsample,mean = para[1],sd=para[2])
  pp=exp(para[3]+para[4]*phenotype)/(1+exp(para[3]+para[4]*phenotype))
  ss=ifelse(runif(Nsample)-pp<0,1,0)
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


    mutc1.hot <- rsparsematrix(k1, Nsample.ps, nnz=rbinom(1, Nsample.ps * k1, hotpp2), rand.x=NULL)
    mutc1.reg <- rsparsematrix(k2, Nsample.ps, nnz=rbinom(1, Nsample.ps * k2, pp2), rand.x=NULL)
    mutc1=rbind(mutc1.hot,mutc1.reg)

    mutc2.hot <- rsparsematrix(k3, Nsample.ps, nnz=rbinom(1, Nsample.ps * k3, hotpp3), rand.x=NULL)
    mutc2.reg <- rsparsematrix(k4, Nsample.ps, nnz=rbinom(1, Nsample.ps * k4, pp3), rand.x=NULL)
    mutc2=rbind(mutc2.hot,mutc2.reg)

    mut.ps <- rbind(mutc1, mutc2)

    mutn1.hot <- rsparsematrix(k1, Nsample.neu, nnz=rbinom(1, Nsample.neu * k1, hotpp1), rand.x=NULL)
    mutn1.reg <- rsparsematrix(k2, Nsample.neu, nnz=rbinom(1, Nsample.neu * k2, pp1), rand.x=NULL)
    mutn1=rbind(mutn1.hot,mutn1.reg)


    mutn2.hot <- rsparsematrix(k3, Nsample.neu, nnz=rbinom(1, Nsample.neu * k3, hotpp1), rand.x=NULL)
    mutn2.reg <- rsparsematrix(k4, Nsample.neu, nnz=rbinom(1, Nsample.neu * k4, pp1), rand.x=NULL)
    mutn2=rbind(mutn2.hot,mutn2.reg)

    mut.neu <- rbind(mutn1, mutn2)

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
  simdata <- list("mutlist"= mutlist, "hotspots"=hotsize, "pheno" = phenotype, "annodata" = annodata, "bmrpars" = bmrpars, "bmrmtxlist" = bmrmtxlist, "para"=para, "hotsize"=hotsize,"phenotype"=phenotype,  "efsize" = list( "betaf0" = betaf0,  "beta_gc" = beta_gc, "avbetaf1" = avbetaf1, "avbetaf2" = avbetaf2, "avbetaf1f2" = avbetaf1f2,"betaf1f2"=betaf1f2))
  return(simdata)
}

