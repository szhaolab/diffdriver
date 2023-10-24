devtools::load_all("~/SimingLab/jiezhou/diffdriver")
parmName=names(parmASHmean$TSG)
## tsg
tsgParm=parmASHmean$TSG
for (i in parmName) {
  j=ifelse (i=="beta_f0","(Intercept)",i)
x=unique(sgmatrix[[1]][,..j])
tt=table(sgmatrix[[1]][,..j])
if (length(tt)==1){
tsgParm[i]=log(2*exp(x*parmASHmean$TSG[i])-1)/x
}else{
  a1=log(2*exp(x[1]*parmASHmean$TSG[i])-1)/x[1]
  a2=log(2*exp(x[2]*parmASHmean$TSG[i])-1)/x[2]
  w=tt[1]/sum(tt)
  tsgParm[i]=a1*w+a2*(1-w)
}
}
names(tsgParm)[6]="(Intercept)"
usethis::use_data(tsgParm,overwrite = T)
## og
ogParm=parmASHmean$OC
for (i in parmName) {
  j=ifelse (i=="beta_f0","(Intercept)",i)
  x=unique(sgmatrix[[1]][,..j])
  tt=table(sgmatrix[[1]][,..j])
  if (length(tt)==1){
    ogParm[i]=log(2*exp(x*parmASHmean$OC[i])-1)/x
  }else{
    a1=log(2*exp(x[1]*parmASHmean$OC[i])-1)/x[1]
    a2=log(2*exp(x[2]*parmASHmean$OC[i])-1)/x[2]
    w=tt[1]/sum(tt)
    ogParm[i]=a1*w+a2*(1-w)
  }
}
names(ogParm)[6]="(Intercept)"
usethis::use_data(ogParm,overwrite = T)


