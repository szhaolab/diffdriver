#!/usr/bin/env Rscript
#library(diffdriver)
library(logging)
  library(data.table)
   library(resample)
   library(Matrix)
   library(fastTopics)
   library(plyr)
  source("ddmodel.R")
  source("comparators.R")
  source("plots.R")
  source("drivermaps_readdatatest.R")
  source("matrixlistToBMR.R")
  source("diffdriver_reg.R")
  source("diffdriver_sig.R")
  source("diffdriver.R")
  

#drivermapsdir <- "~/SimingLab/library/diffdriver_anno/"
drivermapsdir<- "~/SimingLab/library/diffdriver_anno"
outputdir <- "~/temp/"
outputname <- "diffDriver_demo"
load(paste0(drivermapsdir,"/mapsparameters/Fpars.rda"))  # see data-raw folder for how to generate Fpars from drivermaps results 
load(paste0(drivermapsdir,"/mapsparameters/BMRlist.rda")) # see data-raw folder for how to generate BMRlist from drivermaps results

genef <- system.file("extdata/genes.txt", package = "diffdriver")
mutf <- system.file("extdata/mutations.txt", package = "diffdriver")
phenof <- system.file("extdata/phenotypes.txt", package = "diffdriver")
#genef = paste0(drivermapsdir,"/extdata/genes.txt")
#mutf = paste0(drivermapsdir,"/extdata/mutations.txt")
#phenof =paste0(drivermapsdir, "/extdata/phenotypes.txt")

logfile <- file(paste0(outputdir,"/", outputname, ".log"), open="wt")
addHandler(writeToFile, file= logfile, level='DEBUG')
loginfo("Started running diffdriver ...")

res <- diffdriver(genef, mutf, phenof, drivermapsdir = drivermapsdir,mode=2, outputdir = outputdir, outputname = outputname)

save(res,file="res.Rd")


