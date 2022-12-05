
library("devtools")
library(foreach)
load_all("./")

Nsample=100
set.seed(11)
binary=T
bmrpars=log(BMR)
betaf0=1
beta_gc=c(1,0)
para=c(0.8,0.2)

 ###
 simdata <- simulate_2funcv(binary,sgdata, bmrpars, betaf0, Nsample, beta_gc, para,hotseq,hmm)
 ssgdata=simdata$annodata
 mut <- do.call(rbind, simdata$mutlist)
 bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
 hotsize <- do.call(c,simdata$hotsize)
 e <- simdata$pheno


 funcv <- unlist(lapply(ssgdata, "[[", "functypecode"))
 ef <- simdata$efsize
 fe1 <- c(ef$beta_gc[1], ef$beta_gc[1] + ef$beta_gc[2])[as.factor(funcv)]+hotsize
 fe2<- c(ef$avbetaf1, ef$avbetaf1 + ef$avbetaf2)[as.factor(funcv)]
 fe3 <- rep(ef$betaf1f2, length(funcv))
 fe4 <- rep(ef$avbetaf1f2, length(funcv))
 mr <- bmrmtx + ef$betaf0
 foreach(fe=c(fe1,fe2,fe3,fe4)) %dopar%{
   print(paste0("Iteration:",fe))
   res.ddmodel <-  ddmodel(mut,e, mr, fe)
 }
 res <-  ddmodel(mut,e, mr, fe1)

 ##
 simures=power_compare(binary=TRUE,Niter = 1,sgdata=sgdata,Nsample=100,para=c(0.8,0.2),bmrpars=log(BMR),betaf0=1,beta_gc=c(1,0),hotseq=hotseq0,hmm=hmm)


