#!/usr/bin/env Rscript
library("diffdriver")
library("logging")
drivermapsdir<- "/dartfs/rc/lab/S/Szhao/library/diffdriver_anno/"
outputdir <- "~/temp/"
outputname <- "diffDriver_demo"
load(paste0(drivermapsdir,"/mapsparameters/Fpars.rda"))  # see data-raw folder for how to generate Fpars from drivermaps results
load(paste0(drivermapsdir,"/mapsparameters/BMRlist.rda")) # see data-raw folder for how to generate BMRlist from drivermaps results

#genef <- system.file("extdata/genes.txt", package = "diffdriver")
#mutf <- system.file("extdata/mutations.txt", package = "diffdriver")
#phenof <- system.file("extdata/phenotypes.txt", package = "diffdriver")
#genef = paste0(drivermapsdir,"extdata/genes.txt")
#mutf = paste0(drivermapsdir,"extdata/mutations.txt")
#phenof =paste0(drivermapsdir, "extdata/phenotypes.txt")
type="UCEC"
mutationdir="/dartfs/rc/lab/S/Szhao/diffDriver/data_run_prs/tumor_specific_input/"
genef = paste0(drivermapsdir,"extdata/genes.txt")
mutf = paste0(mutationdir,type,"/",type,"_mutations.txt")
phenof =paste0(mutationdir,type,"/",type,"_PRS_46phenotype.txt")
logfile <- file(paste0(outputdir,"/", outputname, ".log"), open="wt")
addHandler(writeToFile, file= logfile, level='DEBUG')
loginfo("Started running diffdriver ...")

res <- diffdriver(genef, mutf, phenof,j=6, drivermapsdir = drivermapsdir,k=6,mode=2, outputdir = outputdir, outputname = outputname)
#save(res,file="res.Rd")


