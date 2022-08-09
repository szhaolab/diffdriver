N=100
aa=c()
for (j in 1:N) {
  seqt=sample(x=c(0,1),size=1,prob =c(hmm[1],hmm[2]))
  for (i in 1:3473) {
    a=ifelse(seqt[i-1]==0,sample(x=c(0,1),size=1,prob =c(hmm[3],hmm[4])),
             sample(x=c(0,1),size=1,prob =c(hmm[5],hmm[6])))
    seqt=c(seqt,a)
  }
  aa=c(aa,sum(seqt))
}

mean(aa)
sd(aa)
