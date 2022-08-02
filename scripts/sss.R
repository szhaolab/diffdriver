source("comparators.R")
source("ddmodel.R")
source("plots.R")
source("simulate_binary.R")
source("simulate_cont.R")
source("simufun.R")
Nsim=1
ncore=2

data(BMR)
data(Fe)
data(sgdata)
#load(file="C:/Users/Jie Zhou/Documents/paper02052022/siming_lab/diffdriver-main/data/diffDriver_demo.ALK.Rd")
hmm=readRDS("hmmOGpar_ASHmean.rds")
hotspot=hmm
hotspot[9]=0
#cl <- makeCluster(ncore,outfile="")
#registerDoParallel(cl)
#print(paste0("start parallel computing using ", ncore, " cores ..."))
    i1=1
    i2=1.2
    i3=400
    family="binary"
    Niter=Nsim
    sgdata=sgdata
    Nsample=i3
    bmrpars=log(BMR)
    betaf0=i1
    beta_gc=c(i2,Fe)
    hotspot= hotspot
set.seed(10)
#simures <-power_compare(family="binary",Niter=Nsim, sgdata=sgdata, bmrpars=bmrpars, Nsample=Nsample, betaf0=betaf0,  beta_gc=beta_gc, para = c(0.8,0.2),hotspot=hotspot)
    simures <-power_compare(family="cont", Niter=Nsim, sgdata=sgdata, Nsample=400, bmrpars=log(BMR), betaf0=i1, beta_gc=c(0,Fe),para=c(0,5,0.1,1),hotspot= hotspot)
    #colnames(simures)=c("beta_gc[1]","beta_gc[2]","old parameter","new parameter")
    #save(simures, file=paste0("power_betaf0=",i1,"_betagc=",i2, "_sample",i3,".Rd"))

