
simulate_1funcv <- function(binary=F,sganno,sgmatrix, bmrpars, betaf0=2, Nsample, beta_gc, para,bpara=c(0,3),tau,rho,hot=0, hmm){
	if (binary==T){ # generate binary phenotype
		Nsamplec <- round(Nsample/2) # number of samples with phenotype E=1 (the rest will be 0)
		Nsamplen <- Nsample-Nsamplec
		phenotype <- c(rep(1,Nsamplec),rep(0,Nsamplen))
		ss=ifelse(phenotype==1,sample(c(0,1),size=Nsamplec,replace=T,prob = c(1-para[1],para[1])),sample(c(0,1),size=Nsamplen,replace=T,prob = c(1-para[2],para[2])))
    Nsample.ps=sum(ss)
    Nsample.neu=Nsample-Nsample.ps
		index=which(ss==1)
		phenotype=c(phenotype[index],phenotype[-index])# the positive samples are placed in the front.
		#u=rnorm(n=length(phenotype),sd=sqrt((1/rho^2-1)*var(phenotype)))+phenotype
		#u1=u/(sum(u^2))^0.5
		#bmrfold=exp(u1*tau)
		bpp=exp(bpara[1]+bpara[2]*phenotype)/(1+exp(bpara[1]+bpara[2]*phenotype))
		bb=ifelse(runif(Nsample)-bpp<0,1,0) # generate pseudo phenotype samples
		}else{
		phenotype=rnorm(Nsample,mean = para[1],sd=para[2])
		pp=exp(para[3]+para[4]*phenotype)/(1+exp(para[3]+para[4]*phenotype))
		bpp=exp(bpara[1]+bpara[2]*phenotype)/(1+exp(bpara[1]+bpara[2]*phenotype))
		ss=ifelse(runif(Nsample)-pp<0,1,0) # generate positive samples
		bb=ifelse(runif(Nsample)-bpp<0,1,0) # generate pseudo phenotype samples
		Nsamplec <- sum(bb)
		Nsamplen <- Nsample-Nsamplec
		pseudophenotype=c(rep(1,Nsamplec),rep(0,Nsamplen))
		index=which(ss==1)
		phenotype=c(phenotype[index],phenotype[-index])# the positive samples are placed in the front.
		Nsample.ps1=sum(ss[1:Nsamplec]) # number of positive samples
		Nsample.ps0=sum(ss[(Nsamplec+1):Nsample]) # number of positive samples
		Nsample.ps=Nsample.ps0+Nsample.ps1
		Nsample.neu1 <- Nsamplec - Nsample.ps1 # number of neutral samples
		Nsample.neu0 <- Nsamplen - Nsample.ps0 # number of neutral samples
		Nsample.neu = Nsample.neu0+ Nsample.neu1
		pseudophenotype=c(pseudophenotype[index],pseudophenotype[-index])# the positive samples are placed in the front.
}
		hotseq=hotspotseq(hmm,sganno) # The first column 'start' is the positions, and the second column 'seqt',is the hotspot indicator.
		if (hot==0){
			hotseq[,2]=0
		}
browser()
	mutlist <- list() # a list of nine mutation matrices
	countlist <- list()
	annodata <- list() # a list of nine annotation data frames
	bmrmtxlist <- list() # a list of nine background mutation matrices
	betagc=c(beta_gc,hmm[9])
	mutRate <- list()
	foldlist <- list()
	for (t in 1:length(sganno)) {
	  hotseqt=merge(sganno[[t]],hotseq,by="start")$seqt
	  selename=names(beta_gc)
		ssgdata=cbind(sgmatrix[[t]][,..selename],hotseqt)
		hotindex=which(hotseqt==1)
		pp.neu=matrix(rep(1,nrow(ssgdata)),ncol=1)%x%matrix(exp(bmrpars[t])*exp(betaf0)*bmrfold,nrow=1)
		fold=as.vector(exp(as.matrix(ssgdata)%*%betagc))
		fold[hotindex]=exp(hmm[9])
		if (any(2*fold<1)){stop("Error:inappropriate parameter settings!")}
		fold=2*fold-1
		pp.ps=ifelse(diag(fold) %*% pp.neu<1,diag(fold)%*%pp.neu,0.99)
		foldlist[[t]]=data.table(fold=fold)
		mutps=apply(pp.ps[,1:Nsample.ps],2,rbinom,n=nrow(pp.ps),size=1)
		mutneu=apply(pp.neu[,(Nsample.ps+1):Nsample],2,rbinom, n=nrow(pp.neu),size=1)

		if (!is.matrix(mutps)){
		mutlist[[t]]=as(cbind(t(mutps),t(mutneu)),"sparseMatrix")
		}else{
		mutlist[[t]]=as(cbind(mutps,mutneu),"sparseMatrix")
		}
		countlist[[t]] <- c(sum(mutlist[[t]]))
		bmrmtxlist[[t]] <- pp.neu
}
# The forllowings are the ture parameters (???)
	fold <- do.call(rbind,foldlist)
	avFe <- rep(log(mean(fold[[1]])*Nsample.ps/Nsample + Nsample.neu/Nsample),nrow(fold))
	diffFe <-  log(fold[[1]]*Nsample.ps/Nsample + Nsample.neu/Nsample)

	simdata <- list("mutlist"= mutlist, "pheno" = phenotype,"foldlist"=fold,"bmrfold"=bmrfold, "annodata" = sganno, "bmrpars" = bmrpars, "bmrmtxlist" = bmrmtxlist, "para"=para, "efsize" = list( "betaf0" = betaf0,  "beta_gc" = betagc, "avFe" = avFe, "diffFe" = diffFe),"nsample"=c(Nsample.ps,Nsample.neu))
	return(simdata)
}





