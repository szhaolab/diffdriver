library("devtools")
library("diffdriver")
library("foreach")
library("doParallel")
source("~/SimingLab/jiezhou/diffdriver/scripts/simufunother.R")
source("~/SimingLab/jiezhou/diffdriver/scripts/simulate.R")
i1=0
load("~/SimingLab/jiezhou/diffdata/gene239/ri200.Rd")
load("~/SimingLab/jiezhou/diffdata/gene239/fanno200.Rd")
i3=200
beta_gc=unlist(parmASHmean[[1]])
names(beta_gc)[6]="(Intercept)"
par1=0
par2=5
par3=0.1
par4=1
set.seed(1)
ncore=10
simuresother=vector("list",200)
cl <- makeCluster(ncore,outfile="")
registerDoParallel(cl)
simuresother1=foreach (j=1:50,.combine='cbind',.packages = c("Matrix", "data.table","diffdriver")) %dopar% {
power_compareother(binary=F,Nite=10,sganno=ri200[[j]],sgmatrix=fanno200[[j]],bmrpars=log(BMR),Nsample=i3,betaf0=i1,beta_gc=beta_gc,para=c(par1,par2,par3,par4),hot=0,hmm=hmm)
}



simuresother2=foreach (j=51:200,.combine='c',.packages = c("Matrix", "data.table","diffdriver")) %dopar% {
power_compareother(binary=F,Nite=100,sganno=ri200[[j]],sgmatrix=fanno200[[j]],bmrpars=log(BMR),Nsample=i3,betaf0=i1,beta_gc=beta_gc,para=c(0,1,0,0),hot=0,hmm=hmm)[[2]]
}
stopCluster(cl)
save(simuresother,file="pppower_betaf0=0_sample200_0.Rd")
q("no")
