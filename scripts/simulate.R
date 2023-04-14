
simulate_1funcv <- function(binary=F,sgdata, bmrpars,faIndex, betaf0=2, Nsample, beta_gc, para,hot=0, hmm){
	if (binary==T){ # generate binary phenotype
		Nsamplec <- round(Nsample/2) # number of samples with phenotype E=1 (the rest will be 0)
		Nsamplen <- Nsample-Nsamplec
		phenotype <- c(rep(1,Nsamplec),rep(0,Nsamplen))
		ss=ifelse(phenotype==1,sample(c(0,1),size=Nsamplec,replace=T,prob = c(1-para[1],para[1])),sample(c(0,1),size=Nsamplen,replace=T,prob = c(1-para[2],para[2])))
		}else{
		phenotype=rnorm(Nsample,mean = para[1],sd=para[2])
		pp=exp(para[3]+para[4]*phenotype)/(1+exp(para[3]+para[4]*phenotype))
		ss=ifelse(runif(Nsample)-pp<0,1,0) # generate positive samples
		}
		index=which(ss==1)
		phenotype=c(phenotype[index],phenotype[-index])# the positive samples are placed in the front.
		Nsample.ps=sum(ss) # number of positive samples
		Nsample.neu <- Nsample - Nsample.ps # number of neutral samples
		hotseq=hotspotseq(hmm,sgdata) # The first column 'start' is the positions, and the second column 'seqt',is the hotspot indicator.
		if (hot==0){
			hotseq[,2]=0
			}
	mutlist <- list() # a list of nine mutation matrices
	countlist <- list()
	annodata <- list() # a list of nine annotation data frames
	bmrmtxlist <- list() # a list of nine background mutation matrices
	betagc=c(beta_gc,hmm[9])
	mutRate <- list()
	foldlist <- list()
	for (t in 1:length(sgdata)) {
	  hotseqt=merge(fanno,hotseq,by=start)$hotseq
		ssgdata=cbind(sgdata[[t]][,c(..faIndex)],hotseqt)
		pp.neu=rep(exp(bmrpars[t])*exp(betaf0),nrow(ssgdata))
		fold=exp(as.matrix(ssgdata)%*%betagc)
		pp.ps=pp.neu*fold
		foldlist[[t]]=data.table(fold=fold)
		mutps=replicate(Nsample.ps,rbinom(length(pp.ps),size=1,pp.ps))
		mutneu=replicate(Nsample.neu,rbinom(length(pp.neu),size=1,pp.neu))
		mutlist[[t]]=as(cbind(mutps,mutneu),"sparseMatrix")
		countlist[[t]] <- c(sum(mutlist[[t]]))
		bmrmtxlist[[t]] <- matrix(exp(bmrpars[t])*exp(betaf0), ncol = ncol(mutlist[[t]]), nrow = nrow(mutlist[[t]])) # background mutation matrix for nytype=t
}

# The forllowings are the ture parameters (???)
	fold <- do.call(rbind,foldlist)
	avFe <- rep(log(mean(fold)*Nsample.ps/Nsample + Nsample.neu/Nsample),nrow(fold))
	diffFe <-  log(fold*Nsample.ps/Nsample + Nsample.neu/Nsample) 

	simdata <- list("mutlist"= mutlist, "pheno" = phenotype,"foldlist"=fold, "annodata" = sgdata, "bmrpars" = bmrpars, "bmrmtxlist" = bmrmtxlist, "para"=para, "efsize" = list( "betaf0" = betaf0,  "beta_gc" = beta_gc, "avFe" = avFe, "diffFe" = diffFe),"nsample"=c(Nsample.ps,Nsample.neu))
	return(simdata)
}





