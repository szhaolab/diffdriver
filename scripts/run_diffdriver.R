#!/usr/bin/env Rscript
library(diffdriver)
library(logging)
drivermapsdir<- "~/SimingLab/library/diffdriver_anno/"
outputdir <- "~/temp/"
outputname <- "diffDriver_demo"
load(paste0(drivermapsdir,"/mapsparameters/Fpars.rda"))  # see data-raw folder for how to generate Fpars from drivermaps results
load(paste0(drivermapsdir,"/mapsparameters/BMRlist.rda")) # see data-raw folder for how to generate BMRlist from drivermaps results

genef <- system.file("extdata/genes.txt", package = "diffdriver")
mutf <- system.file("extdata/mutations.txt", package = "diffdriver")
phenof <- system.file("extdata/phenotypes.txt", package = "diffdriver")

logfile <- file(paste0(outputdir,"/", outputname, ".log"), open="wt")
addHandler(writeToFile, file= logfile, level='DEBUG')
loginfo("Started running diffdriver ...")

res <- diffdriver(genef, mutf, phenof, drivermapsdir = drivermapsdir,mode=2, outputdir = outputdir, outputname = outputname)

save(res,file="res.Rd")


