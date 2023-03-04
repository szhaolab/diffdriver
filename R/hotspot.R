
#' @param The vector hmm is the probabilities for generating the hotspot sequence x(t), in which
#  p(x(0)=0)=hmm[1],p(x(0)=1)=hmm[2],p(x(t)=0|x(t-1)=0)=hmm[3], p(x(t)=1|x(t-1)=0)=hmm[4],
#  p(x(t)=0|x(t-1)=1)=hmm[5],p(x(t)=1|x(t-1)=1)=hmm[6],
#' @param The argument sgdata is the annotation data which provide the position information 
#' @return A n by 2  data frame. The first column start provide the position information, the second column seqt indicates whether this position being a hotspot.  
#' @export
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
