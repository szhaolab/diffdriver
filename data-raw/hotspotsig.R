
load("~/SimingLab/jiezhou/diffdata/gene239/sgannosig.Rd")
sganno=split(sgannosig,sgannosig$nttypecode)
hotspot1sig=vector("list",96)
hotspot2sig=vector("list",96)
for (i in 1:96){
n=nrow(sganno[[i]])
hotspot1sig[[i]]=rep(0,n)
hotspot2sig[[i]]=rep(0,n)
}


## For TCTC type, add hotspots
hotspot1sig[[69]][20]=1
hotspot2sig[[69]][20]=1

hotspot2sig[[69]][50]=1

usethis::use_data(hotspot1sig,overwrite = T)

usethis::use_data(hotspot2sig,overwrite = T)
