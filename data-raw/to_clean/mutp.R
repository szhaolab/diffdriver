
library(devtools)
load_all("~/SimingLab/jiezhou/diffdriver/")
drivermapsdir <- "~/SimingLab/jiezhou/driverMAPS/data/"
Totalnttype <- 9
#Adirbase < paste0(drivermapsdir, "quicktest_data/")
Adirbase <- drivermapsdir
Afileinfo <- list(file = paste(Adirbase, "nttypeXXX_annodata.txt", sep=""),
                  header = c("chrom","start","end","ref","alt","genename","functypecode","nttypecode","expr","repl","hic","mycons","sift","phylop100","MA","ssp","wggerp"),
                  coltype = c("character","numeric","numeric","character","character","character","character","factor","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))

for (j in 1:Totalnttype){
#  dataj <- ddmread_j(Afileinfo, j, varlist = c("chrom","start","genename","functypecode","nttypecode","expr","repl","hic","mycons","sift","phylop100","MA","ssp","wggerp"))
dataj <- ddmread_j(Afileinfo, j, varlist = c("chrom","start","genename","functypecode"))
nttypecode=BMR[j]
geffect=1
functypecodeo=ifelse(dataj[,functypecode]==7,1,
		ifelse(dataj[,functypecode]==8,1,1))
functypecodet=ifelse(dataj[,functypecode]==7,1,
		ifelse(dataj[,functypecode]==8,1,1))
expr=1
repl=1
hic=1
mycons=1
sift=1
phylop100=1
MA=1
ssp=1
wggerp=1
hotspot=1
dataj=cbind(dataj,functypecodeo,functypecodet,nttypecode,geffect,hotspot,expr,repl,hic,mycons,sift,phylop100,MA,ssp,wggerp)
#write.table(dataj,file=paste0("~/SimingLab/jiezhou/driverMAPS/data/nttype",j,"_annorate.txt"))
}

