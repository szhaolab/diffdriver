power_compareother <- function(binary, Niter, sgdata,faIndex=4, Nsample,para,bmrpars,betaf0,beta_gc,hot=0,hmm){
	if (length(faIndex)!=length(beta_gc)){stop("The number of functional annotions does not match!")}
  m1.pvalue <- m2.pvalue <- m3.pvalue <- m4.pvalue <- rep(1,Niter)
	a=c()
	b=c()
	sgdata=split(annoAll$fannomatrixAll[[gene]],annoAll$fannoAll[[gene]]$nttypecode)
	fanno=split(annoAll$fannoAll[[gene]],annoAll$fannoAll[[gene]]$nttypecode)
	for (iter in 1:Niter) {
		simdata <- simulate_1funcv(binary=binary,sgdata, bmrpars,faIndex, betaf0, Nsample, beta_gc, para,hot,hmm)
		ssgdata=simdata$annodata
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		e <- simdata$pheno
		e_bisect=ifelse(e>mean(e),1,0)
		ef <- simdata$efsize
		if (sum(mut) ==0) {next}
		res.m1 <- mlr(mut,e)
		res.m2 <- genefisher(mut,e_bisect)
		res.m3 <- genebinom(mut,e_bisect)
		res.m4 <- genelr(mut,e_bisect)
		m1.pvalue[iter] <-  res.m1$pvalue
		m2.pvalue[iter] <-  res.m2$pvalue
		m3.pvalue[iter] <-  res.m3$pvalue
		m4.pvalue[iter] <-  res.m4$pvalue
		parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		a=rbind(a,parameters)
		nummut=sum(mut)
		b=c(b,nummut)
		}
		return(list("parameters"=a, "m1.pvalue" =m1.pvalue, "m2.pvalue" =m2.pvalue,
               "m3.pvalue" =m3.pvalue,"m4.pvalue" =m4.pvalue,"#mut"=b))
		}



