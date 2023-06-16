load('~/SimingLab/jiezhou/diffdata/gene239/fannoallg.Rd')
load('~/SimingLab/jiezhou/diffdata/gene239/riallg.Rd')
genename200=read.table("~/SimingLab/jiezhou/diffdriver/data-raw/genename200.txt")

fanno200=fannoallg
ri200=riallg


for (i in 1:200) {
mmatrix=fannoallg[[i]]
manno=riallg[[i]][,-c(1,2,3)]
 index1=which(is.na(manno),arr.ind = T)
 index2=which(is.na(mmatrix),arr.ind = T)
 index=c(index1[,1],index2[,1])
#mmatrix=mmatrix[-index,]
#manno=manno[-index,]
fanno200[[i]]=split(mmatrix,manno$nttypecode)
ri200[[i]]=split(manno,manno$nttypecode)
}

save(fanno200,file="~/SimingLab/jiezhou/diffdata/gene239/fanno200.Rd")
save(ri200,file="~/SimingLab/jiezhou/diffdata/gene239/ri200.Rd")
#usethis::use_data(fanno200,overwrite=T)
#usethis::use_data(ri200,overwrite=T)
