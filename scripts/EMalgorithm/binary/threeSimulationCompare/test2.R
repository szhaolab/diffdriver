library(devtools)
load_all("../../../..")
## parameters setting
set.seed(119)
binary=T
bmrpars=log(BMR)
betaf0=1
beta_gc=c(1,0)
para=c(0.8,0.2)
N=100

rr3=c()
rr3_approx=c()
 for (i in 1:N) {
   simdata1 <- simulate_3funcv(sgdata, log(BMR), betaf0 , Nsample = 1000, beta_gc, fracc=0.8, fracn=0.2)
   mut1 <- do.call(rbind, simdata1$mutlist)
   bmrmtx1 <- do.call(rbind, simdata1$bmrmtxlist)
   e1 <- simdata1$pheno
   funcv1 <- unlist(lapply(simdata1$annodata, "[[", "functypecode"))
   ef1 <- simdata1$efsize
   mr1 <- bmrmtx1 + ef$betaf0
   fe1_1 <- c(ef1$avbetaf1,  ef1$avbetaf1+ef1$avbetaf2)[as.factor(funcv1)]
   fe2_1 <- rep(ef1$avbetaf1f2, length(funcv1))
   rr3=c(rr3,ddmodel(mut1,e1, mr1, fe1_1)$pvalue)
rr3_approx=c(rr3_approx,ddmodel_binary_simple(mut1,e1,mr1,fe1_1)$pvalue)
}
summary(rr3)

save(rr3,file="rr3.Rd")
save(rr3_approx,file="rr3_approx.Rd")
