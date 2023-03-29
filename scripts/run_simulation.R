#!/usr/bin/env Rscript
library('devtools')
load_all("../")
library("logging")
drivermapsdir<- "/dartfs/rc/lab/S/Szhao/library/diffdriver_anno/"
outputdir <- "~/temp/"
outputname <- "diffDriver_demo"
type="BLCA"
mutationdir="/dartfs/rc/lab/S/Szhao/diffDriver/data_run_prs/tumor_specific_input/"

genef = paste0(mutationdir,type,"/",type,"genes.txt")
mutf = paste0(mutationdir,type,"/",type,"_mutations.txt")
phenof =paste0(mutationdir,type,"/",type,"_PRS_46phenotype.txt")
hotf=paste0(drivermapsdir,"hotspot.txt")


logfile <- file(paste0(outputdir, outputname, ".log"), open="wt")
addHandler(writeToFile, file= logfile, level='DEBUG')
loginfo("Started running diffdriver ...")

res <- annoAllg(genef, mutf, phenof,j=6,hotf, drivermapsdir = drivermapsdir,k=6,mode=1, outputdir = outputdir, outputname = outputname)
	save(res,file="res.Rd")


