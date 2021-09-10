#!/usr/bin/env Rscript
library(diffdriver)
drivermapsdir <- "~/cancer_somatic/maps/"
load_drivermaps_func(drivermapsdir)

args = commandArgs(trailingOnly=TRUE) # arg is the runconfig file, e.g. data_run/TCGA_2020-01-01/runconfig_TCGA_2020-01-01.R
source(args[1])

log <- file(paste0(args[1],".log"), open="wt")
sink(log)
sink(log, type="message")
print(paste0("Started running diffDriver at ",Sys.time()))

res <- diffDriver(Genef, Mutf, Phenof)

meth <- c("dd", "mlr", "mlr.v2", "fisher", "binom", "lr")
resdf <- data.frame(sapply(meth, function(x)unlist(lapply(lapply(res,'[[',x),'[[','pvalue'))))
colnames(resdf) <- paste0(meth,".p")
resdf[ , paste0(meth,".fdr")] <- apply(resdf,2,p.adjust, method = "fdr")
resdf[,c("mut.E1", "mut.E0", "E1", "E0")] <- do.call(rbind, lapply(lapply(res, '[[', "fisher"), '[[',"count"))

save(res, file = paste0(Outputbase,".Rd"))
write.table(resdf, file= paste0(Outputbase,".txt") , row.names=T, col.names=T, sep="\t", quote = F)
warnings()
print(paste0("Finished at ", Sys.time()))




