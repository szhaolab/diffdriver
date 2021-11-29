# print("simulate mutation in one gene using given parameters, using one functional annotation: LoF. Allow different fractions of samples under selection in E=1 and E=0 groups.")
# print(Sys.time())
library(foreach)
library(doParallel)
library(diffdriver)

#' bmrpars need to be in log scale
power_compare <- function(Niter=200, sgdata, bmrpars, ...){
  m1.pvalue <- m2.pvalue <- m3.pvalue <- m4.pvalue <-
  m5.pvalue <- m6.pvalue <- rep(1,Niter)

  for (iter in 1:Niter) {
    print(paste0("Iteration: ",  iter))
    simdata <- simulate_1funcv(sgdata, bmrpars, ...)
    mut <- do.call(rbind, simdata$mutlist)
    bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
    e <- simdata$pheno
    if (sum(mut) ==0) {next}
    res.m1 <- mlr(mut,e)
    res.m2 <- genefisher(mut,e)
    res.m3 <- genebinom(mut,e)
    res.m4 <- genelr(mut,e)
    funcv <- unlist(lapply(sgdata, "[[", "functypecode"))
    ef <- simdata$efsize
    fe1 <- c(ef$beta_gc[1], ef$beta_gc[1] + ef$beta_gc[2])[as.factor(funcv)]
    mr <- bmrmtx + ef$betaf0
    res.m5 <-  ddmodel(mut,e, mr, fe1)
    fe2 <- rep(ef$avbetaf1f2, length(funcv))
    res.m6 <-  ddmodel(mut,e, mr, fe2)
    m1.pvalue[iter] <-  res.m1$pvalue
    m2.pvalue[iter] <-  res.m2$pvalue
    m3.pvalue[iter] <-  res.m3$pvalue
    m4.pvalue[iter] <-  res.m4$pvalue
    m5.pvalue[iter] <-  res.m5$pvalue
    m6.pvalue[iter] <-  res.m6$pvalue
  }
  return(list( "m1.pvalue" =m1.pvalue, "m2.pvalue" =m2.pvalue,"m3.pvalue" =m3.pvalue,"m4.pvalue" =m4.pvalue,
               "m5.pvalue" =m5.pvalue,"m6.pvalue" =m6.pvalue))
}

Nsim=500
ncore=12

data(BMR)
data(Fe)
data(sgdata)

cl <- makeCluster(ncore,outfile="")
registerDoParallel(cl)
print(paste0("start parallel computing using ", ncore, " cores ..."))

foreach(i1=c(0, 1),.packages = c("Matrix", "data.table")) %:%
  foreach(i2=c(0, 1.2),.packages = "Matrix") %:%
    foreach(i3=c(400,800,1600),.packages = "Matrix")  %dopar% {
      print(c(i1,i2,i3))
      simures <- power_compare(Nsim, sgdata, log(BMR), betaf0=i1, Nsample=i3, beta_gc=c(i2,Fe), fracc=0.8, fracn=0.2)
      save(simures, file=paste0("power_betaf0=",i1,"_betagc=",i2, "_sample",i3,".Rd"))
    }
print("end parallel computing...")
stopCluster(cl)

## plot results
mylist <- list()
for (i1 in c(0, 1)){
  for (i2 in c(0, 1.2)){
    #png(paste0("power_betaf0=",i1,"_betagc=",i2,".png"), 500, 350, units="px", res=300)
    setEPS()
    postscript(file=paste0("power_betaf0=",i1,"_betagc=",i2,".eps"), width=6, height=4.5)
    n<-0
    for(i3 in c(400, 800, 1600)){
      n <- n+1
      load(paste0("power_betaf0=",i1,"_betagc=",i2, "_sample",i3,".Rd"))
      m2.pvalue <- simures[["m2.pvalue"]]
      m4.pvalue <- simures[["m4.pvalue"]]
      m6.pvalue <- simures[["m5.pvalue"]]
      m7.pvalue <- simures[["m6.pvalue"]]
      mylist[[n]] <- c(length(m2.pvalue[m2.pvalue <0.01]),length(m4.pvalue[m4.pvalue <0.01]), length(m6.pvalue[m6.pvalue <0.01]), length(m7.pvalue[m7.pvalue <0.01]))
    }
    outdt <- do.call(cbind, mylist)
    colnames(outdt) <- c(400, 800, 1600)
    mycol=c("#2b8cbe","#fdae61","#FF3D2E","darkgreen")
    plot(1:dim(outdt)[2],ylim=c(0,max(outdt)*1.1), col="white", ylab = "Power", xlab="Number of samples", yaxt="n", xaxt="n")
    yticks_val <- pretty(outdt)
    axis(2, at=yticks_val, lab=paste0(yticks_val/Nsim * 100, " %"),las=2)
    axis(1, at=1:dim(outdt)[2], lab=colnames(outdt))
    grid()
    for (r in 1:4){
      points(outdt[r,],col=mycol[r], pch=19, cex=1.3)
      lines(outdt[r,],col=mycol[r])
    }
    legend("topleft", legend = c("LogisticR","Fisher","diffDriver","diffDriver-baseline"), col= mycol, lwd=1, bty="n", lty=rep(NA,4), pch=rep(19,4))
    dev.off()
  }
}



