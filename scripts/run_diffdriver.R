#!/usr/bin/env Rscript
library(diffdriver)
# You need first download and install drivermaps and run drivermaps for your data before running diffdriver
Drivermapsdir <- "~/cancer_somatic/maps/"

Outputdir <- "~/temp/"
Outputname <- "diffDriver_demo"
data(Fpars)  # see data-raw folder for how to generate Fpars from drivermaps results
data(BMRlist) # see data-raw folder for how to generate BMRlist from drivermaps results

Genef <- system.file("extdata/genes.txt", package = "diffdriver")
Mutf <- system.file("extdata/mutations.txt", package = "diffdriver")
Phenof <- system.file("extdata/phenotypes.txt", package = "diffdriver")

log <- file(paste0(Outputdir,"/", Outputname, ".log"), open="wt")
sink(log)
sink(log, type="message")
print(paste0("Started running diffDriver at ",Sys.time()))

res <- diffdriver(Genef, Mutf, Phenof, drivermapsdir = Drivermapsdir, outputdir = Outputdir, outputname = Outputname)

meth <- c("dd", "mlr", "mlr.v2", "fisher", "binom", "lr")
resdf <- data.frame(sapply(meth, function(x)unlist(lapply(lapply(res,'[[',x),'[[','pvalue'))))
colnames(resdf) <- paste0(meth,".p")
resdf[ , paste0(meth,".fdr")] <- apply(resdf,2,p.adjust, method = "fdr")
resdf[,c("mut.E1", "mut.E0", "E1", "E0")] <- do.call(rbind, lapply(lapply(res, '[[', "fisher"), '[[',"count"))

save(res, file = paste0(Outputdir, "/", Outputname,,".Rd"))
write.table(resdf, file= paste0(Outputdir,"/", Outputname, ".txt") , row.names=T, col.names=T, sep="\t", quote = F)
warnings()
print(paste0("Finished at ", Sys.time()))




