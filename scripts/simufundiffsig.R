
power_comparediff <- function(binary, Niter, sganno,sgmatrix,Nsample,para,signatures,rho,beta_gc,hot=0,hmm,sc){
  m1.pvalue <- rep(1,Niter)
	a=c()
	b=c()

	for (iter in 1:Niter) {
	  print(iter)
		simdata <- simulate_1funcv(binary=binary,sganno=sganno,sgmatrix=sgmatrix, Nsample=Nsample, beta_gc=beta_gc, para=para,signatures=signatures, rho,hot=hot,hmm=hmm,sc=sc)
		ssgdata=simdata$annodata
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		e <- simdata$pheno
		ef <- simdata$efsize
		fe <- ef$diffFe
		mr <- bmrmtx
		if (sum(mut) ==0) {next}
		if (binary == F){
		res.m1 <- ddmodel(mut,e, mr, fe)
		}else{
		res.m1 <- ddmodel_binary_simple(mut,e,mr,fe)
		}
		m1.pvalue[iter] <-  res.m1$pvalue
		parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		a=rbind(a,parameters)
		nummut=sum(mut)
		b=c(b,nummut)
  }
   return(list("parameters"=a, "m1.pvalue" =m1.pvalue,"#mut"=b))
  }



