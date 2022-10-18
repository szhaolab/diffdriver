

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
matrixlistToBMR=function(adirbase,mutf,BMRlist,k=5){
annolist=vector("list",96)
for (i in 1:96) {
  print(paste0("read annotation file: ",i))
annoi=fread(paste0(adirbase,"/TCGA-UCS_nttype",i,"_annodata.txt"))
#write.table(annoi,file=paste0("anno_ACT_C_",i,".txt"))
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
fit1 =fastTopics::fit_poisson_nmf(ymatrix,k = k,numiter = 100,verbose = "none")
#quickplot(x = 1:100,y = fit1$progress$loglik,geom = c("point","line"),
#         xlab = "iteration",ylab = "loglik") + theme_cowplot(12)

fit.multinom <- fastTopics::poisson2multinom(fit1)
## smoothing
#sm=fit.multinom$L%*%t(fit1$F)
ll=fit1$L
colnames(ll)=paste("weight",1:k, sep = "")
ff=fit1$F
colnames(ff)=paste("factor",1:k,sep="")
##
sigmapping=read.table(paste0(adirbase,"config_annotation.txt"))[,-2]
sigmapping[1:10,]
v11=strsplit(sigmapping[,1],split = ",")
muttype=data.table()
for (i in 1:length(v11)) {
  a=v11[[i]]
  aa=paste(substr(a[[1]],1,1),a[[2]],substr(a[[1]],3,3),sep = "")
  bb=paste(a[[2]],">",substr(a[[1]],2,2), sep = "")
  cc=data.table(Type=bb,Subtype=aa)
  muttype=rbind(muttype,cc)
}
sigmapping=cbind(muttype,sigmapping[,2])
sigmapping=sigmapping[substr(Subtype,2,2)=="C" | substr(Subtype,2,2)=="T",]
sigmapping=sigmapping[order(V2),]
colnames(sigmapping)=c("context","alt_allele","index")
ff=cbind(sigmapping[,c("context","alt_allele")],ff)
sm=ll%*%t(ff[,-c(1,2)])
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
ccovariate=plyr::join(covariate,subgenebmr)
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
return(list(sampelsig=wsm,normconstant=normconstant,lambda=genebmr,ll=ll,ff=ff))
}












