

#' Title
#'
#' @param adirbase The path to annotation files
#' @param mutf The path to the mutation files
#' @param BMRlist The parameters from drivermaps
#' @param target The target positions where the bmr's are computed
#'
#' @return
#' @export
#'
#' @examples
matrixlistToBMR=function(adirbase,mutf,BMRlist){
annolist=vector("list",96)
for (i in 1:96) {
  print(paste0("read annotation file: ",i))
annoi=fread(paste0(adirbase,"/TCGA-UCS_nttype",i,"_annodata.txt"))
#annolist[[i]]=annoi[functypecode==6 & ssp==0]
annolist[[i]]=annoi[functypecode==6]
rm(annoi)
}
anno=do.call(rbind,annolist)
rm(annolist)
mutfile <- fread(mutf, header = T)
id= unique(mutfile$SampleID)
nsample=length(id)
ymatrix=matrix(0,nrow = nsample, ncol = 96)
for (i in 1:nsample) {
  idata=mutfile[SampleID==id[i]]
  ni=nrow(idata)
  for (ii in 1:ni) {
    Position=idata$Position[ii]
    Ref=idata$Ref[ii]
    Alt=idata$Alt[ii]
    Chromosome=paste("chr",idata$Chromosome[ii],sep="")
    index=anno[start==Position & ref==Ref &alt==Alt & chrom==Chromosome ,]$nttypecode
    if (length(index)!=1){
      ll=length(index)
      print(paste("sample:",i,"mutation:",ii,"hits:", ll, sep = " "))
    }else{
      ymatrix[i,index]=ymatrix[i,index]+1
    }
  }
}





#the proportion of silent mutation
#sum(ymatrix)/nrow(mutfile)
#save(ymatrix,file = "signature.Rd")
fit1 = fit_poisson_nmf(ymatrix,k = 5,numiter = 100,verbose = "none")
#quickplot(x = 1:100,y = fit1$progress$loglik,geom = c("point","line"),
#         xlab = "iteration",ylab = "loglik") + theme_cowplot(12)

fit.multinom <- poisson2multinom(fit1)
## smoothing
#sm=fit.multinom$L%*%t(fit1$F)
sm=fit1$L%*%t(fit1$F)
## compute the number of silent mutation for each of 96 mutation types.
N=c()
for (i in 1:96) {
  N[i]=anno[functypecode==6 & nttypecode==i,.N,]
}
index=which(N==0)
W=N
W[index]=1
## weighted smoothing
wsm=matrix(0,nrow=nsample,ncol = 96)
wsm=sm%*%diag(1/W)
##gene level
genebmr=data.table(genename=BMRlist$Y_g_s_all$agg_var,lambda=(BMRlist$Y_g_s_all$y+
                                                                    BMRlist$BMpars$fullpars[14])/(BMRlist$Mu_g_s_all$V1+
                                                                                                        BMRlist$BMpars$fullpars[14]))
##position level bmr
coe=BMRlist$BMpars$fullpars[c("expr","repl","hic")]
alpha=BMRlist$BMpars$fullpars[c("alpha")]
covariate=anno[,.(chrom,start,genename,nttypecode,expr,repl,hic)]
rm(anno)
covariate=covariate[expr!="." & repl!="." & hic!=".",]
covariate$expr=as.numeric(covariate$expr)
covariate$repl=as.numeric(covariate$repl)
covariate$hic=as.numeric(covariate$hic)
smooth=covariate$expr*coe[1]+covariate$repl*coe[2]+covariate$hic*coe[3]
covariate=cbind(covariate,linec=smooth)
##
genes=unique(covariate$genename)
lambda=c()
for (i in 1:length(genes)) {
  index=which(genebmr$genename==genes[i])
  aa=ifelse(length(index)!=0,genebmr$lambda[index],1)
  lambda=c(lambda,aa)
}
subgenebmr=data.table(genename= genes,lambda)
ccovariate=join(covariate,subgenebmr)
normconstant=mean(exp(ccovariate$linec)*ccovariate$lambda)
# ee=cbind(ee,explinec=exp(ee$linec))
# ee=cbind(ee,weight=(ee$explinec)*(ee$lambda)/mean((ee$explinec)*ee$lambda))
#rm(anno)
##sample level
# bmrsig=data.table()
# for (i in 1:nsample) {
#   print(paste0("bmr for sampel",i,""))
#   aa=wsm[i,ee$nttypecode]
#   subject=log(ee$weight*aa)
#   bmrsig=cbind(bmrsig,subject)
# }
return(list(sampelsig=wsm,normconstant=normconstant,lambda=genebmr))
}












