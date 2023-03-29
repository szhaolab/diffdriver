
power_comparediff <- function(binary, Niter, sgdata,faIndex=4, Nsample,para,bmrpars,betaf0,beta_gc,hot=0,hmm){
	m1.pvalue <- rep(1,Niter)
	a=c()
	b=c()
	for (iter in 1:Niter) {
		simdata <- simulate_1funcv(binary=binary,sgdata, bmrpars, faIndex,betaf0, Nsample, beta_gc, para,hot,hmm)
		ssgdata=simdata$annodata
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		e <- simdata$pheno
		funcv <- unlist(lapply(ssgdata, "[[", "functypecode"))
		ef <- simdata$efsize
		fe <- c(ef$avbetaf1, ef$avbetaf1 + ef$avbetaf2)[as.factor(funcv)]
		mr <- bmrmtx + ef$betaf0
		if (sum(mut) ==0) {next}
		#res.m1 <- ddmodel(mut,e, mr, fe[[1]])
		res.m1 <- ddmodel_binary_simple(mut,e,mr,fe)
		m1.pvalue[iter] <-  res.m1$pvalue
		parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		a=rbind(a,parameters)
		nummut=sum(mut)
		b=c(b,nummut)
  }
   return(list("parameters"=a, "m1.pvalue" =m1.pvalue,"#mut"=b))
  }



