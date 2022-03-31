# print("simulate mutation in one gene using given parameters, using one functional annotation: LoF. Allow different fractions of samples under selection in E=1 and E=0 groups.")
# print(Sys.time())
library(Matrix)
library(data.table)
library(foreach)
library(doParallel)
codedir <- "/home/simingz/cancer_pheno/cancer_pheno/"
source(paste0(codedir,"comparator_methods.R"))
source(paste0(codedir,"countLRT.R"))

power_compare <- function(sgdata,betaf0=0.5,Nsample=1000,Nc=200,beta_gc=c(0,1), fracc=0.8, fracn=0.2){
  # betaf0, shift of mutation rate from BMR, log scale, shared in all samples.
  # Nsample, total number of samples
  # Nc number of positive genes (associated with phenotype)
  # beta_gc, effect size for positives, log scale
  # effect size for neutral is always 0.
  # fracc, fraction of positively selected samples in group E=1 
  # fracn, fraction of positively selected samples in group E=0
  Nsamplec <- round(Nsample/2) # number of samples with phenotype E=1 (the rest will be 0)
  Nsamplen <- Nsample-Nsamplec
  edata <- c(rep(1,Nsamplec),rep(0,Nsamplen))
  
  Nsamplec.ps <- rbinom(1, Nsamplec, fracc) # number of positively selected samples in group E=1
  Nsamplen.ps <- rbinom(1, Nsamplen, fracn) # number of positively selected samples in group E=0
  Nsample.ps <- Nsamplec.ps + Nsamplen.ps
  Nsample.neu <- Nsample - Nsample.ps
  
  avbetaf0 <- log(exp(betaf0) * (exp(beta_gc[1]) * Nsample.ps/Nsample + Nsample.neu/Nsample)) # average betaf0 for dataset that will be plugged in for cmodel.frac
  avbetaf1 <- log(exp(betaf0) * (exp(beta_gc[1] + beta_gc[2]) * Nsample.ps/Nsample + Nsample.neu/Nsample)) - avbetaf0  # average effect sizes for dataset that will be plugged in for cmodel.frac
  
  m1.pvalue <- rep(1,Nc)
  m2.pvalue <- rep(1,Nc)
  m3.pvalue <- rep(1,Nc)
  m4.pvalue <- rep(1,Nc)
  m5.pvalue <- rep(1,Nc)
  m6.pvalue <- rep(1,Nc)
  m7.pvalue <- rep(1,Nc)
  for (iterg in 1:Nc) {
    mdlist <- list()
    countlist <- list()
    for (t in 1:length(sgdata)){ # Simulate mutation data. t: nucleotide change type
      tnpos1 <- dim(sgdata[[t]][functypecode==7])[1]
      tnpos2 <- dim(sgdata[[t]][functypecode==8])[1]
      mutc1 <- rsparsematrix(tnpos1, Nsample.ps, nnz=rbinom(1, Nsample.ps * tnpos1,BMR[t]*exp(betaf0)*exp(beta_gc[1])), rand.x=NULL)
      mutc2 <- rsparsematrix(tnpos2, Nsample.ps, nnz=rbinom(1, Nsample.ps * tnpos2,BMR[t]*exp(betaf0)*exp(beta_gc[1] + beta_gc[2])), rand.x=NULL)
      mutc <- rbind(mutc1, mutc2)
      mutn1 <- rsparsematrix(tnpos1, Nsample.neu, nnz=rbinom(1, Nsample.neu * tnpos1,BMR[t]*exp(betaf0)), rand.x=NULL)
      mutn2 <- rsparsematrix(tnpos2, Nsample.neu, nnz=rbinom(1, Nsample.neu * tnpos2,BMR[t]*exp(betaf0)), rand.x=NULL)
      mutn <- rbind(mutn1, mutn2)
  
      mutc.out <- cbind(mutc[, 1:Nsamplec.ps], mutn[, 1:(Nsamplec-Nsamplec.ps)])
      mutn.out <- cbind(mutc[, (Nsamplec.ps + 1):Nsample.ps], mutn[, (Nsamplec-Nsamplec.ps+1):Nsample.neu])
      
      mdlist[[t]] <- cbind(mutc.out,mutn.out)
      countlist[[t]] <- c(tnpos1, tnpos2,sum(mutc1),sum(mutc2),sum(mutn1), sum(mutn2), sum(mutc.out), sum(mutn.out), sum(mutc.out[1:tnpos1,]), sum(mutn.out[1:tnpos1,]))
    }
    gmut <- do.call(rbind,mdlist)
    pos1pos2ratio <- colSums(do.call(rbind, countlist))[1]/colSums(do.call(rbind, countlist))[2]
    avbetaf0f1 <- avbetaf0 + log((pos1pos2ratio + exp(avbetaf1))/(pos1pos2ratio+1))
    # print(cmodel.frac(mdlist, edata, sgdata, c(avbetaf0,avbetaf1)))
    # print(cmodel(mdlist, edata))}
    if (sum(gmut) ==0) {next}
    res.m1 <- mlr(mdlist,edata)
    res.m2 <- genefisher(mdlist,edata)
    res.m3 <- genebinom(mdlist,edata)
    res.m4 <- genelr(mdlist,edata)
    res.m5 <- cmodel(mdlist,edata)
    res.m6 <- cmodel.frac(mdlist, edata, sgdata, c(avbetaf0, avbetaf1))
    res.m7 <- cmodel.frac(mdlist, edata, sgdata, c(avbetaf0f1, 0))
    m1.pvalue[iterg] <-  res.m1$coefficients[2,4]
    m2.pvalue[iterg] <-  res.m2
    m3.pvalue[iterg] <-  res.m3
    m4.pvalue[iterg] <-  res.m4$coefficients[2,4]
    m5.pvalue[iterg] <-  res.m5$pvalue
    m6.pvalue[iterg] <-  res.m6$pvalue
    m7.pvalue[iterg] <-  res.m7$pvalue
  }
  return(list( "m1.pvalue" =m1.pvalue, "m2.pvalue" =m2.pvalue,"m3.pvalue" =m3.pvalue,"m4.pvalue" =m4.pvalue,
               "m5.pvalue" =m5.pvalue,"m6.pvalue" =m6.pvalue,"m7.pvalue" =m7.pvalue))
}

BMparsfile <- paste0("~/cancer_somatic/data_run/combined_20170526_5/UCS","/","UCS","_parameters_BMvar.Rdata")
load(BMparsfile)
Totalnttype <- 9
BMR <- exp(BMpars$fullpars[1:Totalnttype])/50
paramdir <- "/project/mstephens/cancer_somatic/maps/param/"
load(paste0(paramdir, "parmASHmean.Rdata")); load(paste0(paramdir, "colmu_sd_funct78.Rdata"))
Fe <- parmASHmean$TSG["functypecode8"]/allsd["functypecode8"] # effect sizes
Fe <- log(exp(Fe)*2)

# readin for one single gene
sg <- "ERBB3"
source("~/cancer_somatic/cancer_somatic/code/R00_config_func.R")
Adirbase <-("~/cancer_somatic/maps/quicktest_data/")
Afileinfo <- list(file = paste(Adirbase, "nttypeXXX_annodata.txt", sep=""),
                  header = c("chrom","start","end","ref","alt","genename","functypecode","nttypecode","expr","repl","hic","mycons","sift","phylop100","MA","ssp","wggerp"),
                  coltype = c("character","numeric","numeric","character","character","character","character","factor","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))

dataall <- list()
sgdata <- list()
for (j in 1:Totalnttype){
  dataall[[j]] <- ddmread_j(Afileinfo, j, varlist = c("chrom","start","genename","functypecode","nttypecode"))
  sgdata[[j]] <- dataall[[j]][(functypecode==7 | functypecode==8)& genename == sg]
}

# save(sgdata, file=paste0("/home/simingz/cancer_pheno/data_run/simulation_2019-05-06/sgdata.Rd"))
# load(paste0("/home/simingz/cancer_pheno/data_run/simulation_2019-04-29/sgdata.Rd"))

Nsim=500

ncore=12
cl <- makeCluster(ncore,outfile="")
registerDoParallel(cl)
print(paste0("start parallel computing using ", ncore, " cores ..."))

foreach(i1=c(0, 1),.packages = c("Matrix", "data.table")) %:%
  foreach(i2=c(0, 1.2),.packages = "Matrix") %:%
    foreach(i3=c(400,800,1600),.packages = "Matrix")  %dopar% {
      print(c(i1,i2,i3))
      simures <- power_compare(sgdata,betaf0=i1,Nsample=i3,Nc=Nsim, beta_gc=c(i2,Fe), fracc=0.8, fracn=0.2)
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
    for(i3 in c(300,500,700,900,1200,1500,2000,3000,4000)){
      n <- n+1
      load(paste0("power_betaf0=",i1,"_betagc=",i2, "_sample",i3,".Rd"))
      m2.pvalue <- simures[["m2.pvalue"]]
      m4.pvalue <- simures[["m4.pvalue"]]
      m6.pvalue <- simures[["m6.pvalue"]]
      m7.pvalue <- simures[["m7.pvalue"]]
      mylist[[n]] <- c(length(m2.pvalue[m2.pvalue <0.01]),length(m4.pvalue[m4.pvalue <0.01]), length(m6.pvalue[m6.pvalue <0.01]), length(m7.pvalue[m7.pvalue <0.01]))
    }
    outdt <- do.call(cbind, mylist)
    colnames(outdt) <- c(300,500,700,900,1200,1500,2000,3000,4000)
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



