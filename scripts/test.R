library("diffdriver")
 i1=0
 i2=1.2
 i3=400
 Nsim=50
 simures=power_compare(family="cont",Niter = Nsim,sgdata=sgdata,bmrpars=log(BMR),Nsample=i3,betaf0=i1,beta_gc=c(i2,Fe),para=c(0,5,0.1,1),hotseq=hotseq0,hmm=hmm)
 save(simures,file="power_betaf0=0_betagc=1.2_sample400.Rd")


 Nsample=i3
 betaf0=i1
