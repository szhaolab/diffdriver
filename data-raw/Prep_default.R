
# TSGpars and OGpars: the default effect sizes for functional annotations.
# This is derived from driverMAPS run results averaging 20 tumor types.


aheader <- c("chrom","start","end","ref","alt","genename","functypecode","nttypecode","expr","repl","hic","mycons","sift","phylop100","MA","ssp","wggerp")

acoltype <- c("character","numeric","numeric","character","character","character","character","factor","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric")

qnvars = c("expr","repl","hic")
bmvars <- c("nttypecode", "expr", "repl", "hic")
bmmuttype <- "functypecode == 6"
funcvars <- c("functypecode", "mycons", "sift", "phylop100", "MA")
functypecodelevel <- "7"
funcvmuttype <- "functypecode == 7 | functypecode == 8"
readinvars <- c("chrom", "start", "ref", "alt", "genename", "ssp", bmvars, funcvars)
bmreadinvars <- c("chrom", "start", "ref", "alt", "genename", "ssp", bmvars, "functypecode")
 # will be normalized, except for nttypecode

load( "inst/extdata/parmASHmean.Rdata")
TSGpars <- parmASHmean[[1]]
OGpars <- parmASHmean[[2]]

OGs <- read.table("inst/extdata/OG.txt",header=F)
TSGs <- read.table("inst/extdata/TSG.txt",header=F)

hmm=readRDS("inst/extdata/hmmOGpar_ASHmean.rds")


# nttype code (96)
sigmapping=read.table("/dartfs-hpc/rc/home/m/f0052zm/szhao_lab_share/library/diffdriver_anno/coding_data_03212024.txt", header = T)
sigmapping <- cbind(sigmapping, do.call(rbind, strsplit(sigmapping[,1],split = ",")))
colnames(sigmapping)[2:4] <- c("index", "context", "alt")
sigmapping$ref <- substring(sigmapping[,3],2,2)
sigmapping <- sigmapping[order(sigmapping$index),]
sigmapping <- sigmapping[seq(1,192,2), c("index", "context", "ref","alt")]


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
                  funcvars, functypecodelevel, funcvmuttype, bmreadinvars,
                  readinvars, qnvars, sigmapping,
                  internal = TRUE,overwrite = TRUE)



