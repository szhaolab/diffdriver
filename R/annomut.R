
#'For each gene and each subject, extract the annotation and mutation data from raw data
#' @param sgdata A list with length 9. Each of its components is a data frame with annotation information. The first eight columns are
#' chromosome,start,end,ref,alt,genename,functypecode,nttypecode respectively. The other columns
#' include other position level information.
#' @param mutation Data frame with mutation information with columns being
#' Chromosome, Position, Ref, Alt, SampleID respectively.
#' @return A list named annomut in which length(annomut)=length(genename). Each compoent of annomut
#' is a list with length equal to the sample size.
bmrate=function(sgdata,mutation){
  ssgdata=do.call(rbind,sgdata)
p1=ncol(ssgdata)
genename=unique(ssgdata$genename)
annomut=vector("list",length = length(genename))
sampleid=unique(mutation$SampleID)
n=length(sampleid)
for (k in 1:length(genename)) {
  annomut[[k]]=vector("list",length = n)
  for (i in 1:n) {
index=which(ssgdata$genename==genename[k])
chromosome=ssgdata$chrom[index[1]]
posit=unique(ssgdata$start[index])
anno1k=data.frame(chromosome,genename[k],posit)
anno2k=data.frame(matrix(0,nrow=nrow(kdata),ncol=p1-6))
colnames(anno2k)=colnames(ssgdata)[-c(1:6)]
annogenek=cbind(anno1k,anno2k)
annomut[[k]][[i]]=rbind(annogene,annogenek)
indexi=which(mutation$SampleID==sampleid[i])
mutationi=mutation[indexi,]
for (ii in 1:nrow(mutationi)) {
  aa=mutation[ii,2:4]
  index1=which(annomut[[k]][[i]]$posit==aa[1])
  index2=which(ssgdata[,3:5]==aa)
  annomut[[k]][[i]][index1,-c(1:3)]=ssgdata[index2,6:p1]
}
}
}
return(annomut)
}
