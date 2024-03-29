# readin data for one single gene

library(devtools)
load_all("~/SimingLab/jiezhou/diffdriver/")
browser()
#drivermapsdir <- "~/cancer_somatic/maps/"
#drivermapsdir <- "~/SimingLab/library/diffdriver_anno/"
drivermapsdir <- "~/SimingLab/jiezhou/driverMAPS/data/"
Totalnttype <- 9
#Adirbase < paste0(drivermapsdir, "quicktest_data/")
Adirbase <- drivermapsdir
Afileinfo <- list(file = paste(Adirbase, "nttypeXXX_annodata.txt", sep=""),
                  header = c("chrom","start","end","ref","alt","genename","functypecode","nttypecode","expr","repl","hic","mycons","sift","phylop100","MA","ssp","wggerp"),
                  coltype = c("character","numeric","numeric","character","character","character","character","factor","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))

dataall <- list()
sgdata78 <- list()
for (j in 1:Totalnttype){
  dataall[[j]] <- ddmread_j(Afileinfo, j, varlist = c("chrom","start","genename","functypecode","nttypecode"))
  sgdata78[[j]] <- dataall[[j]][(functypecode==7 | functypecode==8)]
}
genename200=unique(sgdata78[[1]]$genename)[1:250]
index1=grep(',',genename200)
index2=grep(';',genename200)
index=c(index1,index2)
genename200=genename200[-index]
sgdata200 <- vector("list",length(genename200))
for (i in 1:length(genename200)){
sgdata200[[i]]=vector("list",9)
for (j in 1:Totalnttype){
sgdata200[[i]][[j]]=sgdata78[[j]][genename==genename200[i]]
}
}
names(sgdata200)=genename200
save(sgdata200,file="sgdata200.Rd")
write(genename200,file="genename200.txt")
usethis::use_data(genename200,overwrite= T)
usethis::use_data(sgdata200, overwrite = T)
load('./data-raw/signatures.Rd')
colnames(signatures)[1]="type"
usethis::use_data(signatures,overwrite = T)

load('./data-raw/codeSignature.Rd')

usethis::use_data(codeSignature,overwrite = T)

og=fread("OG.txt",header=F)
usethis::use_data(og,overwrite=T)
