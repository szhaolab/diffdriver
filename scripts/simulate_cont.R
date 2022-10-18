simulate_2funcv <- function(sgdata, bmrpars, betaf0=2, Nsample=1000, beta_gc=c(0.1,0.5), par=c(0,1,0.5,0.5)){
  phenotype=rnorm(Nsample,mean = par[1],sd=par[2])
  pp=exp(par[3]+par[4]*phenotype)/(1+exp(par[3]+par[4]*phenotype))
  ss=ifelse(runif(Nsample)-pp<0,1,0)
  Nsample.ps=sum(ss)
  Nsample.neu <- Nsample - Nsample.ps

  mutlist <- list()
  countlist <- list()
  annodata <- list()
  bmrmtxlist <- list()
  for (t in 1:length(sgdata)) {
    tnpos1 <- dim(sgdata[[t]][functypecode==7])[1]
    tnpos2 <- dim(sgdata[[t]][functypecode==8])[1]
    annodata[[t]] <- rbind(sgdata[[t]][functypecode==7], sgdata[[t]][functypecode==8])
    mut1=matrix(nrow = tnpos1,ncol = Nsample)
    mut2=matrix(nrow = tnpos2,ncol = Nsample)
    pp1=exp(bmrpars[t])*exp(betaf0)*exp(beta_gc[1])
    pp2=exp(bmrpars[t])*exp(betaf0)*exp(beta_gc[1] + beta_gc[2])
    pp3=exp(bmrpars[t])*exp(betaf0)
    for (j in 1:Nsample) {
      if (ss[j]==1){
        mut1[,j]= sample(c(0,1),tnpos1, replace=T,prob= c(1-pp1,pp1))
        mut2[,j]= sample(c(0,1),tnpos2, replace = T, prob = c(1-pp2,pp2))
      }
      if (ss[j]==0){
        mut1[,j]= sample(c(0,1),tnpos1, replace=T, prob= c(1-pp3,pp3))
        mut2[,j]= sample(c(0,1),tnpos2, replace = T, prob = c(1-pp3,pp3))
      }
    }
    mut=rbind(mut1,mut2)
    mutlist[[t]] <- mut
    countlist[[t]] <- c(tnpos1, tnpos2,sum(mut),sum(mut1), sum(mut2))
    bmrmtxlist[[t]] <- matrix(bmrpars[t], ncol = ncol(mutlist[[t]]), nrow = nrow(mutlist[[t]]))
  }

  avbetaf1 <- log(exp(beta_gc[1]) * Nsample.ps/Nsample + Nsample.neu/Nsample)
  avbetaf2 <- log(exp(beta_gc[1] + beta_gc[2]) * Nsample.ps/Nsample + Nsample.neu/Nsample)
  pos1pos2ratio <- colSums(do.call(rbind, countlist))[1]/colSums(do.call(rbind, countlist))[2]
  avbetaf1f2 <- log((pos1pos2ratio*exp(avbetaf1) + exp(avbetaf2))/(pos1pos2ratio+1))
  betaf1f2 <- log((pos1pos2ratio*exp(beta_gc[1]) + exp(beta_gc[2]))/(pos1pos2ratio+1))
  simdata <- list("mutlist"= mutlist, "pheno" = phenotype, "annodata" = annodata, "bmrpars" = bmrpars, "bmrmtxlist" = bmrmtxlist, "par"=par, "efsize" = list( "betaf0" = betaf0,  "beta_gc" = beta_gc, "avbetaf1" = avbetaf1, "avbetaf2" = avbetaf2, "avbetaf1f2" = avbetaf1f2,"betaf1f2"=betaf1f2))
  return(simdata)
}

