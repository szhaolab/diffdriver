library(data.table)
og=fread("OG.txt",header=F)
save(og,file="og.Rd")
usethis::use_data(og,overwrite=T)
