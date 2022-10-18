## This code compare the mutation number for hotspots and non-hotspots
## so that you can get the idea how hotspot cans increase the risk

library(diffdriver)
bmrpars=log(BMR)
betaf0=1
Nsample=10000
beta_gc=c(1.2,Fe)
para1=c(0.8,0.2)
para2=c(0,5,0.1,1)


matsgdata=do.call(rbind,sgdata)
ssgdata=merge(matsgdata,hotseq,by="start")
simdata1 <- simulate_1funcv(sgdata, bmrpars, betaf0, Nsample, beta_gc, para1,hotseq,hmm)
simdata2 <- simulate_2funcv(sgdata, bmrpars, betaf0, Nsample, beta_gc, para2,hotseq,hmm)
hot1=simdata1$hotspots
hot2=simdata2$hotspots
mut1=simdata1$mutlist
mut2=simdata2$mutlist
result1=matrix(nrow = 9,ncol=2)
for (t in 1:length(sgdata)){
index1=which(hot1[[t]]!=0)
hotmut1=mut1[[t]][index1,]
if (length(index1)>0){
  coldmut1=mut1[[t]][-index1,]
  }else{
    coldmut1=mut1[[t]]
}
result1[t,1]=ifelse(length(index1)>0,sum(hotmut1)/length(index1),0)
result1[t,2]=sum(coldmut1)/(nrow(mut1[[t]])-length(index1))
}

result2=matrix(nrow = 9,ncol=2)
for (t in 1:length(sgdata)){
  index2=which(hot2[[t]]!=0)
  hotmut2=mut2[[t]][index2,]
  if (length(index2)>0){
    coldmut2=mut2[[t]][-index2,]
  }else{
    coldmut2=mut2[[t]]
  }
  result2[t,1]=ifelse(length(index2)>0,sum(hotmut2)/length(index2),0)
  result2[t,2]=sum(coldmut2)/(nrow(mut2[[t]])-length(index2))
}

colSums(result1)/Nsample*100
colSums(result2)/Nsample*100

