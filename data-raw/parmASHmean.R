load('parmASHmean.Rdata')
og=read.table("OG1.txt")
tsg=read.table("TSG1.txt")
usethis::use_data(parmASHmean,overwrite=T)
