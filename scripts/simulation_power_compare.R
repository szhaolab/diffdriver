family="cont"
Niter=100
sgdata=sgdata
bmrpars=log(BMR)
betaf0=1
Nsample=400
beta_gc=c(1.2,Fe)
para=c(0,5,0.1,1)
hotseq = hotseq
hmm = hmm
library(Matrix)
library(diffdriver)
library(data.table)
hmm[9]=0
for (i1 in c(0, 1)){
  for (i2 in c(0, 1.2)){
  for (i3 in c(400,800,1600)) {
    print(c(i1,i2,i3))
    simures=power_compare(family="binary", Niter=100, sgdata=sgdata, bmrpars=log(BMR), betaf0=i1,Nsample=i3,beta_gc=c(i2,Fe),para=c(0.8,0.2),hotseq = hotseq,hmm = hmm)
    save(simures,file=paste0("power_betaf0=",i1,"_betagc=",i2,"_sample",i3,".Rd")
)
    }
  }
}

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
      m3.pvalue <- simures[["m3.pvalue"]]
      m4.pvalue <- simures[["m4.pvalue"]]
      mylist[[n]] <- c(length(m1.pvalue[m1.pvalue <0.01]),length(m2.pvalue[m2.pvalue <0.01]),length(m3.pvalue[m3.pvalue <0.01]), length(m4.pvalue[m4.pvalue <0.01]))
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
    for (r in 1:4){
      points(outdt[r,],col=mycol[r], pch=19, cex=1.3)
      lines(outdt[r,],col=mycol[r])
    }
    legend("topleft", legend = c("LinearR","bionomial","LogisticR","Fisher"), col= mycol, lwd=1, bty="n", lty=rep(NA,4), pch=rep(19,4))
    dev.off()
  }
}


