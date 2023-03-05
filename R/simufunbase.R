#' Title
#'
#' @param family indicates the type of E variable
#' @param Niter is the number of replication in simulation
#' @param sgdata is the annotation data
#' @param bmrpars is the background mutation rate.
#' @param ... other parameters
#' @param Nsample, total number of samples
#' @param para para is the parameters for generating the data
#' @param bmrpars a vector, each item is bmr in log scale for one nttype
#' @param betaf0, shift of mutation rate from BMR, log scale, shared in all samples.
#' @param beta_gc, effect size for functional covariates, log scale. Right now we only have one functional covariate, that is whether the mutation is missense or loss of function mutation. beta_gc[1] indicates the shift of mutation rate for missense ( coded as 7 in functype code column in `sgdata`), beta_gc[2] indicates the shift of mutation rate for loss of function ( coded as 8 in functype code column in `sgdata`). log scale.
#' @param hotspot hotspot[1] is the probability of being hotspot for a given position. hotspot[2] is the log size
#' for hotspots.
#' @return A list composed of the p-values for 8 models
#' and the parameters used in these models.
#' @export
power_comparebase <- function(binary, Niter, sgdata, Nsample,para,bmrpars,betaf0,beta_gc,hot=0,hmm){
	m1.pvalue <-  rep(1,Niter)
	a=c()
	b=c()
	for (iter in 1:Niter) {
		simdata <- simulate_1funcv(binary=binary,sgdata, bmrpars, betaf0, Nsample, beta_gc, para,hot,hmm)
		ssgdata=simdata$annodata
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		hotsize <- do.call(c,simdata$hotsize)
		e <- simdata$pheno
		funcv <- unlist(lapply(ssgdata, "[[", "functypecode"))
		ef <- simdata$efsize
		fe <- rep(ef$avbetaf1f2, length(funcv))
		mr <- bmrmtx + ef$betaf0
		if (sum(mut) ==0) {next}
		#res.m1 <- ddmodel(mut,e, mr, fe)
		res.m1 <- ddmodel_binary_simple(mut,e,mr,fe)
		m1.pvalue[iter] <-  res.m1$pvalue
		parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		a=rbind(a,parameters)
		nummut=sum(mut)
		b=c(b,nummut)
  }
	return(list("parameters"=a, "m1.pvalue" =m1.pvalue,"#mut"=b))
  }
