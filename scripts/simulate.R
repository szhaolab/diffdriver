
simulate_1funcv <- function(binary=F,sgdata, bmrpars,faIndex=4, betaf0=2, Nsample, beta_gc, para,hot=0, hmm){
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
	foldMatrix <- list()
	for (t in 1:length(sgdata)) {
		ssgdata=merge(sgdata[[t]][,c(2,..faIndex)],hotseq,by="start")
		ssgdata$functypecode=ifelse(ssgdata$functypecode==7,0,1)
		pp.neu=rep(exp(bmrpars[t])*exp(betaf0),nrow(ssgdata))
		fold=exp(as.matrix(ssgdata[,-1])%*%betagc[-1]+betagc[1])
		pp.ps=pp.neu*fold
		foldMatrix[[t]]=setDT(data.frame(foldNumber=fold,functypecode=ssgdata[,functypecode]))
		tnpos1 <- dim(ssgdata[,functypecode==7])[1]
		tnpos2 <- dim(ssgdata[,functypecode==8])[1]
		mutps=replicate(Nsample.ps,rbinom(length(pp.ps),size=1,pp.ps))
		mutneu=replicate(Nsample.neu,rbinom(length(pp.neu),size=1,pp.neu))
		mutframe=cbind(mutps,mutneu)
		coordinate=which(mutframe==1,arr.ind = T)
		mutlist[[t]]=sparseMatrix(i=coordinate[,1],j=coordinate[,2],dims = c(nrow(ssgdata),Nsample))
		countlist[[t]] <- c(tnpos1, tnpos2,sum(mutlist[[t]]))
		bmrmtxlist[[t]] <- matrix(bmrpars[t], ncol = ncol(mutlist[[t]]), nrow = nrow(mutlist[[t]])) # background mutation matrix for nytype=t
}

# The forllowings are the ture parameters (???)
#	avbetaf1 <- log(exp(beta_gc[1]) * Nsample.ps/Nsample + Nsample.neu/Nsample)
#	avbetaf2 <- log(exp(beta_gc[1] + beta_gc[2]) * Nsample.ps/Nsample + Nsample.neu/Nsample)
#	pos1pos2ratio <- colSums(do.call(rbind, countlist))[1]/colSums(do.call(rbind, countlist))[2]
#	avbetaf1f2 <- log((pos1pos2ratio*exp(avbetaf1) + exp(avbetaf2))/(pos1pos2ratio+1))
	rate <- do.call(rbind,foldMatrix)
	avbetaf1f2 <- log(mean(rate[,foldNumber])*Nsample.ps/Nsample + Nsample.neu/Nsample)
	rate7 <- do.call(rbind,foldMatrix)[functypecode==0,]
	rate8 <- do.call(rbind,foldMatrix)[functypecode==1,]
	avbetaf1 <- log(mean(rate7[,foldNumber]) * Nsample.ps/Nsample + Nsample.neu/Nsample)
	avbetaf2 <- log(mean(rate8[,foldNumber]) * Nsample.ps/Nsample + Nsample.neu/Nsample)

	simdata <- list("mutlist"= mutlist, "pheno" = phenotype,"foldMatrix"=foldMatrix, "annodata" = sgdata, "bmrpars" = bmrpars, "bmrmtxlist" = bmrmtxlist, "para"=para, "efsize" = list( "betaf0" = betaf0,  "beta_gc" = beta_gc, "avbetaf1" = avbetaf1, "avbetaf2" = avbetaf2, "avbetaf1f2" = avbetaf1f2),"nsample"=c(Nsample.ps,Nsample.neu))
	return(simdata)
}





