
# TSGpars and OGpars: the default effect sizes for functional annotations.
# This is derived from driverMAPS run results averaging 20 tumor types.


aheader <- c("chrom","start","end","ref","alt","genename","functypecode","nttypecode","expr","repl","hic","mycons","sift","phylop100","MA","ssp","wggerp")

acoltype <- c("character","numeric","numeric","character","character","character","character","factor","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric")

bmvars <- c("nttypecode", "expr", "repl", "hic")
bmmuttype <- "functypecode == 6"
funcvars <- c("functypecode", "mycons", "sift", "phylop100", "MA")
functypecodelevel <- "7"
funcvmuttype <- "functypecode == 7 | functypecode == 8"
readinvars <- c("genename", "ssp", bmvars, funcvars)
bmreadinvars <- c("genename", "ssp","nttypecode",bmvars, funcvars)
qnvars = c("expr","repl","hic") # will be normalized, except for nttypecode


load( "inst/extdata/parmASHmean.Rdata")
TSGpars <- parmASHmean[[1]]
OGpars <- parmASHmean[[2]]

OGs <- read.table("inst/extdata/OG.txt",header=F)
TSGs <- read.table("inst/extdata/TSG.txt",header=F)

hmm=readRDS("inst/extdata/hmmOGpar_ASHmean.rds")


# nttype code (96)
annodir <- "~/temp/annodir96"
sigmapping=read.table(file.path(annodir,"config_annotation.txt"))[,-2]
sigmapping[1:10,]
v11=strsplit(sigmapping[,1],split = ",")
muttype=data.table()
for (i in 1:length(v11)) {
  a=v11[[i]]
  aa=paste(substr(a[[1]],1,1),a[[2]],substr(a[[1]],3,3),sep = "")
  bb=paste(a[[2]],">",substr(a[[1]],2,2), sep = "")
  cc=data.table(Type=bb,Subtype=aa)
  muttype=rbind(muttype,cc)
}
sigmapping=cbind(muttype,sigmapping[,2])
sigmapping=sigmapping[substr(Subtype,2,2)=="C" | substr(Subtype,2,2)=="T",]
sigmapping=sigmapping[order(V2),]
colnames(sigmapping)=c("context","alt_allele","index")

# ## compute the number of silent mutation for each of 96 mutation types.
# annolist=vector("list", 96)
# for (i in 1:96) {
#   print(paste0("read annotation file: ",i))
#   annoi=fread(paste0(annodir,"/TCGA-UCS_nttype",i,"_annodata.txt"))
#   annolist[[i]]=annoi[functypecode==6]
#   rm(annoi)
# }
# anno=do.call(rbind,annolist)
# rm(annolist)
# W=c()
# for (i in 1:96) {
#   W[i]=anno[functypecode==6 & nttypecode==i,.N,]
# }
# W[which(W==0)] = max(W)

usethis::use_data(OGs, TSGs, OGpars, TSGpars, hmm,
                  aheader, acoltype, bmvars, bmmuttype,
                  funcvars, functypecodelevel, funcvmuttype,
                  readinvars, bmreadinvars, qnvars,
                  sigmapping,
                  internal = TRUE,overwrite = TRUE)



