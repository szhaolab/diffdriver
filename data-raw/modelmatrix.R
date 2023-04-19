load('./data-raw/fannoallg.Rd')
load('./data-raw/riallg.Rd')
#usethis::use_data(fannoallg,overwrite=T)
#usethis::use_data(riallg,overwrite=T)

erbb3matrix=fannoallg$ERBB3
erbb3anno=riallg$ERBB3[,-c(1,2,3)]
 index1=which(is.na(erbb3anno),arr.ind = T)
 index2=which(is.na(erbb3matrix),arr.ind = T)
 index=c(index1[,1],index2[,1])
sganno=split(erbb3anno,erbb3anno$nttypecode)
sgmatrix=split(erbb3matrix,erbb3anno$nttypecode)
usethis::use_data(sganno,overwrite=T)
usethis::use_data(sgmatrix,overwrite=T)
