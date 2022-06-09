#!/usr/bin/env Rscript

#library(diffdriver)
library(logging)
library(data.table)
library(resample)
library(Matrix)
library(fastTopics)
library(plyr)
library(diffdriver)

Drivermapsdir <- "~/cancer_somatic/maps/"
Outputdir <- "~/temp/"
Outputname <- "diffDriver_demo"
data(Fpars)  # see data-raw folder for how to generate Fpars from drivermaps results
data(BMRlist) # see data-raw folder for how to generate BMRlist from drivermaps results

Genef <- system.file("extdata/genes.txt", package = "diffdriver")
Mutf <- system.file("extdata/mutations.txt", package = "diffdriver")
Phenof <- system.file("extdata/phenotypes.txt", package = "diffdriver")

logfile <- file(paste0(Outputdir,"/", Outputname, ".log"), open="wt")
addHandler(writeToFile, file= logfile, level='DEBUG')
loginfo("Started running diffdriver ...")
res <- diffdriver(genef, mutf, phenof, drivermapsdir = drivermapsdir, outputdir = outputdir, outputname = outputname)

# meth <- c("dd", "mlr", "mlr.v2", "fisher", "binom", "lr")
# resdf <- data.frame(sapply(meth, function(x)unlist(lapply(lapply(res,'[[',x),'[[','pvalue'))))
# colnames(resdf) <- paste0(meth,".p")
# resdf[ , paste0(meth,".fdr")] <- apply(resdf,2,p.adjust, method = "fdr")
# resdf[,c("mut.E1", "mut.E0", "E1", "E0")] <- do.call(rbind, lapply(lapply(res, '[[', "fisher"), '[[',"count"))
#
# save(res, file = paste0(Outputdir, "/", Outputname,".Rd"))
# write.table(resdf, file= paste0(Outputdir,"/", Outputname, ".txt") , row.names=T, col.names=T, sep="\t", quote = F)
# warnings()
# loginfo("diffdriver finished.")




