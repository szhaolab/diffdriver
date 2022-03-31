# print("simulate mutation in one gene using given parameters, using one functional annotation: LoF")
# print(Sys.time())
library(Matrix)
library(data.table)
library(foreach)
library(doParallel)
codedir <- "/home/simingz/cancer_pheno/cancer_pheno/"
source(paste0(codedir,"comparator_methods.R"))
source(paste0(codedir,"countLRT.R"))

power_compare <- function(sgdata,betaf0=0.5,Nsample=1000,Nc=200,beta_gc=c(1,1),beta_gn=c(0,0)){
  # Nsample, total number of samples
  # Nc number of positive genes (associated with phenotype)
  # beta_gc, effect size for positives
  # beta_gn, effect size for negatives
  Nsamplec <- round(Nsample/2) # number of samples with phenotype E=1 (the rest will be 0)
  Nsamplen <- Nsample-Nsamplec
  edata <- c(rep(1,Nsamplec),rep(0,Nsamplen))
  
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
      mutc1 <- rsparsematrix(tnpos1, Nsamplec, nnz=rbinom(1,Nsamplec*tnpos1,BMR[t]*exp(betaf0)*exp(beta_gc[1])), rand.x=NULL)
      mutc2 <- rsparsematrix(tnpos2, Nsamplec, nnz=rbinom(1,Nsamplec*tnpos2,BMR[t]*exp(betaf0)*exp(beta_gc[1] + beta_gc[2])), rand.x=NULL)
      mutc <- rbind(mutc1, mutc2)
      mutn1 <- rsparsematrix(tnpos1, Nsamplen, nnz=rbinom(1,Nsamplen*tnpos1,BMR[t]*exp(betaf0)*exp(beta_gn[1])), rand.x=NULL)
      mutn2 <- rsparsematrix(tnpos2, Nsamplen, nnz=rbinom(1,Nsamplen*tnpos2,BMR[t]*exp(betaf0)*exp(beta_gn[1] + beta_gn[2])), rand.x=NULL)
      mutn <- rbind(mutn1, mutn2)
      mdlist[[t]] <- cbind(mutc,mutn)
      countlist[[t]] <- c(tnpos1, tnpos2,sum(mutc1),sum(mutc2),sum(mutn1), sum(mutn2))
    }
    gmut <- do.call(rbind,mdlist)
    print(colSums(do.call(rbind, countlist)))
    if (sum(gmut) ==0) {next}
    res.m1 <- mlr(mdlist,edata)
    res.m2 <- genefisher(mdlist,edata)
    res.m3 <- genebinom(mdlist,edata)
    res.m4 <- genelr(mdlist,edata)
    res.m5 <- cmodel(mdlist,edata)
    res.m6 <- cmodel.fe1(mdlist,edata, Fe, sgdata)
    res.m7 <- cmodel.fe2(mdlist,edata, Fe, sgdata)
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
  sgdata[[j]] <- dataall[[j]][(functypecode==7 | functypecode==8 )& genename == sg]
}

save(sgdata, file=paste0("/home/simingz/cancer_pheno/data_run/simulation_2019-04-29/sgdata.Rd"))
# load(paste0("/home/simingz/cancer_pheno/data_run/simulation_2019-04-29/sgdata.Rd"))

Nsim=5000
# # for loop version
# for (i1 in c(0,1)){
#   for (i2 in c(0, 0.2, 0.8)){
#     for (i3 in c(300, 1000, 3000)){
#       print(c(i1,i2,i3))
#       simures <- power_compare(sgdata,betaf0=i1,Nsample=i3,Nc=Nsim, beta_gc=i2,beta_gn=0)
#       save(simures , file=paste0("power_betaf0=",i1,"_betagc=",i2, "_sample",i3,".Rd"))
#     }
#   }
# } # with Nc=1000, each power_compare run takes around 8min. with small Nc like 100, it is <1 minute.

# foreach loop version
ncore=12
cl <- makeCluster(ncore,outfile="")
registerDoParallel(cl)
print(paste0("start parallel computing using ",ncore, " cores ..."))
foreach(i1=c(0,1),.packages = c("Matrix", "data.table")) %:%
  foreach(i2=c(0, 0.5, 1.2),.packages = "Matrix") %:%
    foreach(i3=c(300, 1000),.packages = "Matrix")  %dopar% {
      print(c(i1,i2,i3))
      if (i2==0) {Ferun <- 0} else {Ferun <- Fe}
      simures <- power_compare(sgdata,betaf0=i1,Nsample=i3,Nc=Nsim, beta_gc=c(i2,Ferun),beta_gn=c(0,0))
      save(simures , file=paste0("power_betaf0=",i1,"_betagc=",i2, "_sample",i3,".Rd"))
    }
print("end parallel computing...")
stopCluster(cl)

## plot results
for (i1 in c(0,1)){
  for (i2 in c(0, 0.5, 1.2)){
    png(paste0("power_betaf0=",i1,"_betagc=",i2,".png"), 1000, 600)
    par(mfrow=c(1,2),mar=c(3,3,3,2.1),oma=c(1,1,5,0))
    for (i3 in c(300, 1000)){
      load(paste0("power_betaf0=",i1,"_betagc=",i2, "_sample",i3,".Rd"))
      m1.pvalue <- simures[["m1.pvalue"]]
      m2.pvalue <- simures[["m2.pvalue"]]
      m3.pvalue <- simures[["m3.pvalue"]]
      m4.pvalue <- simures[["m4.pvalue"]]
      m5.pvalue <- simures[["m5.pvalue"]]
      m6.pvalue <- simures[["m6.pvalue"]]
      m7.pvalue <- simures[["m7.pvalue"]]
      barplot(c(length(m1.pvalue[m1.pvalue <0.01]),length(m2.pvalue[m2.pvalue <0.01]), length(m3.pvalue[m3.pvalue <0.01]),
                length(m4.pvalue[m4.pvalue <0.01]),length(m5.pvalue[m5.pvalue <0.01]),length(m6.pvalue[m6.pvalue <0.01]),
                length(m7.pvalue[m7.pvalue <0.01])), 
                main="Power comparison",col=c("darkgreen","salmon","blue","grey","orange","orangered","orchid"), ylim=c(0,Nsim))
      #axis(side=1, at=seq(0.7,8,1))
      legend("topleft",
             legend = c("ANOVA", "Binomial","LogisticR","Fisher","CountLRT","CountLRT-anno1","CountLRT-anno2"),
             fill = c("darkgreen","salmon","blue","grey","orange","orangered","orchid"))
    }
    mtext('# sample=300',at=.25,side=3,outer=T,cex=1.2)
    mtext('# sample=1000',at=.75,side=3,outer=T,cex=1.2)
    dev.off()
  }
}

