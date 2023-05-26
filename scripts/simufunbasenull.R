
power_comparebasenull <- function(binary, Niter, sganno,sgmatrix, Nsample,para,bmrpars,betaf0,beta_gc,hot=0,hmm){
  m1.pvalue <-  rep(1,Niter)
	a=c()
	b=c()

	for (iter in 1:Niter) {
		simdata <- simulate_1funcv(binary=binary,sganno,sgmatrix, bmrpars, betaf0, Nsample, beta_gc, para,hot,hmm)
		ssgdata=simdata$annodata
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		e <- simdata$pheno
		ef <- simdata$efsize
		fe <- ef$avFe
		mr <- bmrmtx
		momr <- mr*(fe%*%t(rep(1,ncol(mr))))
		if (sum(mut) ==0) {next}
		if (binary==F){
		res.m1 <- ddmodel(mut,e, momr, fe)
		}else{
		res.m1 <- ddmodel_binary_simple(mut,e,momr,fe)
		}
		m1.pvalue[iter] <-  res.m1$pvalue
		parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		a=rbind(a,parameters)
		nummut=sum(mut)
		b=c(b,nummut)
  }
	return(list("parameters"=a, "m1.pvalue" =m1.pvalue,"#mut"=b))
  }




