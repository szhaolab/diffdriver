
simulate_1funcv <- function(binary=F,sganno,sgmatrix, Nsample, beta_gc, para,hot=0, hmm=hmm,signatures,rho,sc){
	if (binary==T){ # generate binary phenotype
		Nsamplec <- round(Nsample/2) # number of samples with phenotype E=1 (the rest will be 0)
		Nsamplen <- Nsample-Nsamplec
		phenotype <- c(rep(1,Nsamplec),rep(0,Nsamplen))
		ss=ifelse(phenotype==1,sample(c(0,1),size=Nsamplec,replace=T,prob = c(1-para[1],para[1])),sample(c(0,1),size=Nsamplen,replace=T,prob = c(1-para[2],para[2])))
		selection=rbind(ss,1-ss)
		Nsample.ps=sum(ss)
   	Nsample.neu=Nsample-Nsample.ps
   	sig1=merge(signatures,codeSignature)
   	sig2=sig1[order(sig1$number),c(3,2,4,5)]
		bmrfold= bmrSignature(e=phenotype,signatures=sig2,rho=rho,sc=sc)# generate the 96 times n bmr matrix
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
		# hotseq=hotspotseq(hmm,sganno) # The first column 'start' is the positions, and the second column 'seqt',is the hotspot indicator.
		# if (hot==0){
		# 	hotseq[,2]=0
		# }

	mutlist <- list() # a list of nine mutation matrices
	countlist <- list()
	annodata <- list() # a list of nine annotation data frames
	bmrmtxlist <- list() # a list of nine background mutation matrices
	betagc=c(beta_gc,hmm[9])
	mutRate <- list()
	foldlist <- list()
	size=c()
	for (t in 1:length(sganno)) {
	  	size=c(size,nrow(sganno[[t]]))
	  	if (nrow(sganno[[t]])==0) {next}
	  	hotseqt= hotspot2sig[[t]]
	  	selename=names(beta_gc)
		ssgdata=cbind(sgmatrix[[t]][,..selename],hotseqt)
		hotindex=which(hotseqt==1)
		pp.neu=matrix(rep(1,nrow(ssgdata)),ncol=1)%x%matrix(bmrfold$bmr[t,],nrow=1)
		fold=as.vector(exp(as.matrix(ssgdata)%*%betagc))
		fold[hotindex]=exp(hmm[9])
		if (any(2*fold<1)){stop("Error:inappropriate parameter settings!")}
		fold=2*fold-1
		F=cbind(fold,1)%*%selection
		pp.total=ifelse(F*pp.neu<1,F*pp.neu,0.99)
		foldlist[[t]]=data.table::data.table(fold=fold)
		if (nrow(pp.total)>1) {
		mutlist[[t]]=as(apply(pp.total,2,rbinom,n=nrow(pp.total),size=1),"sparseMatrix")
		}else{
		  mutlist[[t]]=as(t(rbinom(length(pp.total),size=1,pp.total)),"sparseMatrix")
		  }
		countlist[[t]] <- c(sum(mutlist[[t]]))
		bmrmtxlist[[t]] <- log(pp.neu)
}
# The forllowings are the ture parameters (???)
	fold <- do.call(rbind,foldlist)
	avFe <- rep(log(mean(fold[[1]])*Nsample.ps/Nsample + Nsample.neu/Nsample),nrow(fold))
	diffFe <-  log(fold[[1]]*Nsample.ps/Nsample + Nsample.neu/Nsample)
 	#covariate=apply(bmrfold$bmr, 2, "%*%",size)
	covariate=apply(bmrfold$loadings, 2, sum)
	simdata <- list("mutlist"= mutlist, "pheno" = phenotype,"foldlist"=fold,"covariate"=covariate,"bmrfold"=bmrfold, "annodata" = sganno, "bmrmtxlist" = bmrmtxlist, "para"=para, "efsize" = list(  "avFe" = avFe, "diffFe" = diffFe),"nsample"=c(Nsample.ps,Nsample.neu))
	return(simdata)
}





