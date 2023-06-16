hotspotseq=function(hmm,sgdata){
  pos=sort(unique(do.call(rbind,sgdata)$start))
  seqt=sample(x=c(0,1),size=1,prob =c(hmm[1],hmm[2]))
  for (i in 2:length(pos)) {
    a=ifelse(seqt[i-1]==0,sample(x=c(0,1),size=1,prob =c(hmm[3],hmm[4])),
             sample(x=c(0,1),size=1,prob =c(hmm[5],hmm[6])))
    seqt=c(seqt,a)
  }
  hotseq=cbind(start=pos,seqt=seqt)
  return(hotseq)
}
hotseq=hotspotseq(hmm,sgdata)
sum(hotseq[,2])
hotseq6=hotseq
usethis::use_data(hotseq6, overwrite = T)
load("~/SimingLab/jiezhou/diffdriver/data-raw/nttypeNumber.Rd")
usethis::use_data(N,overwrite = T)
