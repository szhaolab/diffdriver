# print("simulate mutation in one gene using given parameters, using one functional annotation: LoF. Allow different fractions of samples under selection in E=1 and E=0 groups.")
# print(Sys.time())
library(foreach)
library(doParallel)
library(diffdriver)

#' bmrpars need to be in log scale
power_compare <- function(family="binary", Niter=200, sgdata, bmrpars, ...){
  m1.pvalue <- m2.pvalue <- m3.pvalue <- m4.pvalue <-
    m5.pvalue <- m6.pvalue <- rep(1,Niter)
a=c()
  for (iter in 1:Niter) {
    print(paste0("Iteration: ",  iter))
    if (family="binary"){
    simdata <- simulate_1funcv(sgdata, bmrpars, ...)
    }else{
      simdata <- simulate_2funcv(sgdata, bmrpars, ...)
    }
    ssgdata=simdata$annodata
    mut <- do.call(rbind, simdata$mutlist)
    bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
    e <- simdata$pheno
    e_bisect=ifelse(e>median(e),1,0)
    if (sum(mut) ==0) {next}
    res.m1 <- mlr(mut,e)
    res.m2 <- genefisher(mut,e_bisect)
    res.m3 <- genebinom(mut,e_bisect)
    res.m4 <- genelr(mut,e)
    funcv <- unlist(lapply(ssgdata, "[[", "functypecode"))
    ef <- simdata$efsize
    fe1 <- c(ef$beta_gc[1], ef$beta_gc[1] + ef$beta_gc[2])[as.factor(funcv)]
    fe2<- c(ef$avbetaf1, ef$avbetaf1 + ef$avbetaf2)[as.factor(funcv)]
    fe3 <- rep(ef$betaf1f2, length(funcv))
    fe4 <- rep(ef$avbetaf1f2, length(funcv))
    mr <- bmrmtx + ef$betaf0
    res.m5 <-  ddmodel(mut,e, mr, fe1)
    res.m6<- ddmodel(mut,e, mr, fe2)
    parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
    a=rbind(a,parameters)
    res.m7 <-  ddmodel(mut,e, mr, fe3)
    res.m8 <-  ddmodel(mut,e, mr, fe4)
     m1.pvalue[iter] <-  res.m1$pvalue
     m2.pvalue[iter] <-  res.m2$pvalue
     m3.pvalue[iter] <-  res.m3$pvalue
     m4.pvalue[iter] <-  res.m4$pvalue
     m5.pvalue[iter] <-  res.m5$pvalue
     m6.pvalue[iter] <-  res.m6$pvalue
     m7.pvalue[iter] <-  res.m7$pvalue
     m8.pvalue[iter] <-  res.m8$pvalue
  }
  return(list("parameters"=a, "m1.pvalue" =m1.pvalue, "m2.pvalue" =m2.pvalue,"m3.pvalue" =m3.pvalue,"m4.pvalue" =m4.pvalue,
               "m5.pvalue" =m5.pvalue,"m6.pvalue" =m6.pvalue,
              "m7.pvalue" =m7.pvalue,"m8.pvalue" =m8.pvalue))
}

Nsim=100
ncore=2

data(BMR)
data(Fe)
data(sgdata)
load(file="C:/Users/Jie Zhou/Documents/paper02052022/siming_lab/diffdriver-main/data/diffDriver_demo.ALK.Rd")

cl <- makeCluster(ncore,outfile="")
registerDoParallel(cl)
print(paste0("start parallel computing using ", ncore, " cores ..."))
a=c()
foreach(i1=c(0, 1),.packages = c("Matrix", "data.table")) %:%
  foreach(i2=c(0, 1.2),.packages = "Matrix") %:%
  foreach(i3=c(400,800,1600),.packages = c("Matrix", "diffdriver"))  %dopar% {
    print(c(i1,i2,i3))
    simures <-power_compare(Nsim, sgdata, log(BMR), betaf0=i1, Nsample=i3, beta_gc=c(i2,Fe), fracc=0.8, fracn=0.2)
    colnames(simures)=c("beta_gc[1]","beta_gc[2]","old parameter","new parameter")
    #save(simures, file=paste0("power_betaf0=",i1,"_betagc=",i2, "_sample",i3,".Rd"))
    }
stopCluster(cl)

## plot results
mylist <- list()
for (i1 in c(0, 1)){
  for (i2 in c(0, 1.2)){
    png(paste0("power_betaf0=",i1,"_betagc=",i2,".png"), 500, 350,units="px")
    #setEPS()
    #postscript(file=paste0("power_betaf0=",i1,"_betagc=",i2,".eps"), width=6, height=4.5)
    n<-0
    for(i3 in c(400, 800, 1600)){
      n <- n+1
      load(paste0("power_betaf0=",i1,"_betagc=",i2, "_sample",i3,".Rd"))
      m1.pvalue<-  simures[["m1.pvalue"]]
      m2.pvalue <- simures[["m2.pvalue"]]
      m4.pvalue <- simures[["m4.pvalue"]]
      m6.pvalue <- simures[["m5.pvalue"]]
      m7.pvalue <- simures[["m6.pvalue"]]
      mylist[[n]] <- c(length(m1.pvalue[m1.pvalue <0.01]),length(m2.pvalue[m2.pvalue <0.01]),length(m4.pvalue[m4.pvalue <0.01]), length(m6.pvalue[m6.pvalue <0.01]), length(m7.pvalue[m7.pvalue <0.01]))
    }
    outdt <- do.call(cbind, mylist)
    colnames(outdt) <- c(400, 800, 1600)
    mycol=c("#00AFBB","#2b8cbe","#fdae61","#FF3D2E","darkgreen")
    #plot(1:dim(outdt)[2],ylim=c(0,max(outdt)*1.1), col="white", ylab = "Power", xlab="Number of samples", yaxt="n", xaxt="n")
    plot(1:dim(outdt)[2],ylim=c(0,100), col="white", ylab = "Power", xlab="Number of samples", yaxt="n", xaxt="n")
    #yticks_val <- pretty(outdt)
    yticks_val <- seq(0,100,20)
    #axis(2, lab=paste0(yticks_val/Nsim * 100, " %"),las=2)
    axis(2, at=yticks_val, lab=paste0(yticks_val/Nsim * 100, " %"),las=2)
    axis(1, at=1:dim(outdt)[2], lab=colnames(outdt))
    grid()
    for (r in 1:5){
      points(outdt[r,],col=mycol[r], pch=19, cex=1.3)
      lines(outdt[r,],col=mycol[r])
    }
    legend("topleft", legend = c("LinearR","LogisticR","Fisher","diffDriver","diffDriver-baseline"), col= mycol, lwd=1, bty="n", lty=rep(NA,4), pch=rep(19,4))
    dev.off()
  }
}


load(paste0("power_betaf0=",0,"_betagc=",0, "_sample",800,".Rd"))
load("power_betaf0=0_betagc=0_sample400.Rd")
