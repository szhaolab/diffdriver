library("data.table")
library("fastTopics")
library(ggplot2)
library(cowplot)
library(diffdriver)
##read into data
bmrsig=function(gene,m=3){
link1=paste(gene,"/",gene,"_mutations.txt",sep="")
link2=paste(gene,"/BMRlist.Rd",sep = "")
mutlink <- system.file(link1, package = "diffdriver")
paralink= system.file(link2, package = "diffdriver")
ss1link=Paras= system.file("ssgdata1.rda", package = "diffdriver")
ss2link=Paras= system.file("ssgdata2.rda", package = "diffdriver")
mutfile=fread(mutlink)
load(ss1link)
load(ss2link)
load(paralink)
sigpub=as.data.table(read.csv("TCGA_SBS_signatures_patterns.csv"))
sigpub$Subtype=substr(sigpub$Subtype,2,4)
sigmapping=read.table("config_annotation.txt")[,-2]
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
sigpubb=merge(sigmapping,sigpub)
sigpubb=sigpubb[order(sigpubb$V2),]
anno1=do.call(rbind,ssgdata1)[functypecode==6,]
anno2=do.call(rbind,ssgdata2)[functypecode==6,]
rm(ssgdata1)
rm(ssgdata2)
##
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
index1=anno1[start==Position & ref==Ref &alt==Alt & chrom==Chromosome ,]$nttypecode
index2=anno2[start==Position & ref==Ref &alt==Alt & chrom==Chromosome ,]$nttypecode

if ((length(index1) + length(index2))!=1){
  ll=length(index1)+length(index2)
  print(paste("sample:",i,"mutation:",ii,"hits:", ll, sep = " "))
  }else{
  nttypecode=ifelse(length(index1)==1,index1,index2)
  ymatrix[i,nttypecode]=ymatrix[i,nttypecode]+1
}
}
}
#the proportion of silent mutation
#sum(ymatrix)/nrow(mutfile)
#save(ymatrix,file = "signature.Rd")
fit1 = fit_poisson_nmf(ymatrix,k = 5,numiter = 100,verbose = "none")
#quickplot(x = 1:100,y = fit1$progress$loglik,geom = c("point","line"),
  #        xlab = "iteration",ylab = "loglik") + theme_cowplot(12)

fit.multinom <- poisson2multinom(fit1)
## smoothing
sm=fit.multinom$L%*%t(fit1$F)

## compute the number of silent mutation for each of 96 mutation types.
N=c()
for (i in 1:96) {
N[i]=anno1[functypecode==6 & nttypecode==i,.N,]+
  anno2[functypecode==6 & nttypecode==i,.N,]
}
index=which(N==0)
W=N/sum(N)
W[index]=1
## weighted smoothing
wsm=matrix(0,nrow=nsample,ncol = 96)
wsm=sm%*%diag(1/W)
##gene level

genebmr=data.table(genename=BMRlist[[gene]]$Y_g_s_all$agg_var,
        lambda=(BMRlist[[gene]]$Y_g_s_all$y
                +BMRlist[[gene]]$BMpars$fullpars[14])/(BMRlist[[gene]]$Mu_g_s_all$V1+
                                               BMRlist[[gene]]$BMpars$fullpars[14]))
##position level bmr
coe=BMRlist[[gene]]$BMpars$fullpars[c("expr","repl","hic")]
alpha=BMRlist[[gene]]$BMpars$fullpars[c("alpha")]
covariate1=anno1[,.(chrom,start,genename,nttypecode,expr,repl,hic)]
covariate1=covariate1[expr!="." & repl!="." & hic!=".",]
covariate1$expr=as.numeric(covariate1$expr)
covariate1$repl=as.numeric(covariate1$repl)
covariate1$hic=as.numeric(covariate1$hic)
covariate2=anno2[,.(chrom,start,genename, nttypecode,expr,repl,hic)]
covariate2=covariate2[expr!="." & repl!="." & hic!=".",]
covariate2$expr=as.numeric(covariate2$expr)
covariate2$repl=as.numeric(covariate2$repl)
covariate2$hic=as.numeric(covariate2$hic)
smooth1=covariate1$expr*coe[1]+covariate1$repl*coe[2]+covariate1$hic*coe[3]
smooth2=covariate2$expr*coe[1]+covariate2$repl*coe[2]+covariate2$hic*coe[3]
covariate1=cbind(covariate1,linec=smooth1)
covariate2=cbind(covariate2,linec=smooth2)
covariate=rbind(covariate1,covariate2)
rm(anno1,anno2,covariate1,covariate2,smooth1,smooth2)
##
shared=intersect( covariate$genename,genebmr$gene)
covariate=covariate[genename %in% shared,]
genebmr=genebmr[genename %in% shared,]
ee=merge(covariate,genebmr,by="genename")
rm(covariate,genebmr,shared)
ee=cbind(ee,explinec=exp(ee$linec))
ee=cbind(ee,weight=(ee$explinec)*(ee$lambda)/sum(ee$explinec))
##sample level
rm(list=setdiff(ls(),c("wsm","ee","nsample","sigpubb","fit.multinom")))
for (i in 1:nsample) {
  print(i)
  sigbmr=wsm[i,ee$nttypecode]*(ee$weight)
  save(sigbmr,file = paste("sample",i,".Rd",sep = ""))
}
return(list(sigpubb=sigpubb,sigest=fit.multinom$F))
}
result=bmrsig(gene="UCS")
par(mar=c(1,1,1,1))
par(mfrow=c(15,2))
for (i in 4:33) {
plot(1:96,result$sigpubb[,i])
}
par(mfrow=c(5,1))
for (i in 1:5) {
  plot(1:96,result$sigest[,i])
}
