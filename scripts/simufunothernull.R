power_compareother <- function(binary, Niter, sganno,sgmatrix, Nsample,para,bpara=c(0,3),tau=1,rho,bmrpars,betaf0,beta_gc,hot=0,hmm){

  m1.pvalue <- m2.pvalue <- m3.pvalue <- m4.pvalue <- rep(1,Niter)
	a=c()
	b=c()

	for (iter in 1:Niter) {
		simdata <- simulate_1funcv(binary=binary,sganno,sgmatrix, bmrpars, betaf0, Nsample, beta_gc, para,bpara,tau,rho,hot,hmm)
		ssgdata=simdata$annodata
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		bmrfold <- simdata$bmrfold
		e <- simdata$pheno
		e_bisect=ifelse(e>mean(e),1,0)
		ef <- simdata$efsize
		if (sum(mut) ==0) {next}
		res.m1 <- mlr(mut,e)
		res.m2 <- genefisher(mut,e_bisect)
		res.m3 <- genebinom(mut,e_bisect)
		res.m4 <- genelr(mut,e_bisect,covariate= bmrfold)
		m1.pvalue[iter] <-  res.m1$pvalue
		m2.pvalue[iter] <-  res.m2$pvalue
		m3.pvalue[iter] <-  res.m3$pvalue
		m4.pvalue[iter] <-  res.m4$pvalue
		parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		a=rbind(a,parameters)
		nummut=sum(mut)
		b=c(b,nummut)
		}
		return(list("m1.pvalue" =m1.pvalue, "m2.pvalue" =m2.pvalue,
               "m3.pvalue" =m3.pvalue,"m4.pvalue" =m4.pvalue,"#mut"=b))
		}



