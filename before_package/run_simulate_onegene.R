# print("simulate mutation in one gene using given parameters")
# print(Sys.time())
library(Matrix)
library(foreach)
library(doParallel)
source("comparator_methods.R")
source("countLRT.R")

power_compare <- function(sgdata,betaf0=0.5,Nsample=1000,Nc=200,beta_gc=1,beta_gn=0){
  # Nsample, total number of samples
  # Nc number of positive genes (associated with phenotype)
  # beta_gc, effect size for positives
  # beta_gn, effect size for negatives
  Nsamplec <- round(Nsample/2) # number of samples with phenotype E=1 (the rest will be 0)
  Nsamplen <- Nsample-Nsamplec
  edata <- c(rep(1,Nsamplec),rep(0,Nsamplen))

  # m1.beta <- rep(NA,Nc)
  # m2.beta <- matrix(NA,nrow=2,ncol=Nc)
  m1.pvalue <- rep(1,Nc)
  m2.pvalue <- rep(1,Nc)
  m3.pvalue <- rep(1,Nc)
  m4.pvalue <- rep(1,Nc)
  m5.pvalue <- rep(1,Nc)
  for (iterg in 1:Nc) {
    mdlist <- list()
    for (t in 1:length(sgdata)){ # Simulate mutation data. t: nucleotide change type
      tnpos <- dim(sgdata[[t]])[1]
      mutc <- rsparsematrix(tnpos, Nsamplec, nnz=rbinom(1,Nsamplec*tnpos,BMR[t]*exp(betaf0)*exp(beta_gc)), rand.x=NULL)
      mutn <- rsparsematrix(tnpos, Nsamplen, nnz=rbinom(1,Nsamplen*tnpos,BMR[t]*exp(betaf0)*exp(beta_gn)), rand.x=NULL)
      mdlist[[t]] <- cbind(mutc,mutn)
    }
    gmut <- do.call(rbind,mdlist)
    if (sum(gmut) ==0) {next}
    res.m1 <- mlr(mdlist,edata)
    res.m2 <- cmodel(mdlist,edata)
    res.m3 <- genebinom(mdlist,edata)
    res.m4 <- genelr(mdlist,edata)
    res.m5 <- genefisher(mdlist,edata)
    m1.pvalue[iterg] <-  res.m1$coefficients[2,4]
    m2.pvalue[iterg] <-  res.m2$pvalue
    m3.pvalue[iterg] <-  res.m3
    m4.pvalue[iterg] <-  res.m4$coefficients[2,4]
    m5.pvalue[iterg] <-  res.m5
    # m1.beta[iterg] <-  res.m1$coefficients[2,1]
    # m2.beta[,iterg] <-  res.m2$alt.beta
  }
  return(list( "m1.pvalue" =m1.pvalue, "m2.pvalue" =m2.pvalue,"m3.pvalue" =m3.pvalue,"m4.pvalue" =m4.pvalue,"m5.pvalue" =m5.pvalue))
}

BMparsfile <- paste0("~/cancer_somatic/data_run/combined_20170526_5/UCS","/","UCS","_parameters_BMvar.Rdata")
load(BMparsfile)
Totalnttype <- 9
BMR <- exp(BMpars$fullpars[1:Totalnttype])/50

## readin for one single gene
# sg <- "ERBB3"
# source("~/cancer_somatic/cancer_somatic/code/R00_config_func.R")
# Adirbase <-("~/cancer_somatic/maps/sample_data/")
# Afileinfo <- list(Afile = paste(Adirbase, "nttypeXXX_annodata.txt", sep=""),
#                   header = c("chrom","start","end","ref","alt","genename","functypecode","nttypecode","expr","repl","hic","mycons","sift","phylop100","MA","ssp","wggerp"),
#                   coltype = c("character","numeric","numeric","character","character","character","factor","factor","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))
# 
# dataall <- list()
# sgdata <- list()
# for (j in 1:Totalnttype){
#   dataall[[j]] <- ddmread_j(Afileinfo, j, varlist = c("chrom","start","genename","functypecode","nttypecode"))
#   sgdata[[j]] <- dataall[[j]][functypecode==7 & genename == sg]
# }
## or load from a previous run
load("~/cancer_somatic/data_run/simulation_20180628/sgdata.Rd")

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
ncore=18
cl <- makeCluster(ncore,outfile="")
registerDoParallel(cl)
print(paste0("start parallel computing using ",ncore, " cores ..."))
foreach(i1=c(0,1),.packages = "Matrix") %:%
  foreach(i2=c(0, 0.2, 0.8),.packages = "Matrix") %:%
    foreach(i3=c(300, 1000, 3000),.packages = "Matrix")  %dopar% {
      print(c(i1,i2,i3))
      simures <- power_compare(sgdata,betaf0=i1,Nsample=i3,Nc=Nsim, beta_gc=i2,beta_gn=0)
      save(simures , file=paste0("power_betaf0=",i1,"_betagc=",i2, "_sample",i3,".Rd"))
    }
print("end parallel computing...")
stopCluster(cl)

## plot results
for (i1 in c(0,1)){
  for (i2 in c(0, 0.2, 0.8)){
    png(paste0("power_betaf0=",i1,"_betagc=",i2,".png"), 1000, 300)
    par(mfrow=c(1,3),mar=c(3,3,3,2.1),oma=c(1,1,5,0))
    for (i3 in c(300, 1000, 3000)){
      load(paste0("power_betaf0=",i1,"_betagc=",i2, "_sample",i3,".Rd"))
      m1.pvalue <- simures[["m1.pvalue"]]
      m2.pvalue <- simures[["m2.pvalue"]]
      m3.pvalue <- simures[["m3.pvalue"]]
      m4.pvalue <- simures[["m4.pvalue"]]
      m5.pvalue <- simures[["m5.pvalue"]]
      # plot(jitter(rep(1,Nc)),m1.beta,xlim=c(0.9,1.1), ylim=c(-0.2,1), xaxt='n',xlab=NA,ylab=NA,col="darkgreen",
      #      main="Method 1 parameter")
      # axis(side=1, at=1, labels = "effect size")
      # plot(jitter(c(rep(1,Nc),rep(1.3,Nc))),c(m2.beta[1,],m2.beta[2,]),xlim=c(0.9,1.4), ylim=c(-0.5,2),
      #      main="Method 2 parameters", xaxt='n',xlab=NA,ylab=NA,col="salmon")
      # axis(side=1, at=c(1,1.3), labels = c("beta_f0","beta_g"))
      # segments(0.9,betaf0,1.1,betaf0,col="red")
      # segments(1.2,beta_gc,1.4,beta_gc,col="red")
      
      barplot(c(length(m1.pvalue[m1.pvalue <0.01]),length(m2.pvalue[m2.pvalue <0.01]), length(m3.pvalue[m3.pvalue <0.01]),length(m4.pvalue[m4.pvalue <0.01]),length(m5.pvalue[m5.pvalue <0.01])),
              main="Power comparison",col=c("darkgreen","salmon","blue","grey","orange"), ylim=c(0,Nsim))
      axis(side=1, at=seq(0.7,6,1.2))
      legend("topright", 
             legend = c("ANOVA", "CountLRT","Binomial","LogisticR","Fisher"), 
             fill = c("darkgreen","salmon","blue","grey","orange"))
    }
    mtext('# sample=300',at=.16,side=3,outer=T,cex=1.2)
    mtext('# sample=1000',at=.5,side=3,outer=T,cex=1.2)
    mtext('# sample=3000',at=.83,side=3,outer=T,cex=1.2)
    dev.off()
  }
}

