#!/usr/bin/env Rscript
library(diffdriver)
library(logging)
<<<<<<< HEAD
# You need first download and install drivermaps and run drivermaps for your data before running diffdriver
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

res <- diffdriver(Genef, Mutf, Phenof, drivermapsdir = Drivermapsdir, outputdir = Outputdir, outputname = Outputname)

meth <- c("dd", "mlr", "mlr.v2", "fisher", "binom", "lr")
resdf <- data.frame(sapply(meth, function(x)unlist(lapply(lapply(res,'[[',x),'[[','pvalue'))))
colnames(resdf) <- paste0(meth,".p")
resdf[ , paste0(meth,".fdr")] <- apply(resdf,2,p.adjust, method = "fdr")
resdf[,c("mut.E1", "mut.E0", "E1", "E0")] <- do.call(rbind, lapply(lapply(res, '[[', "fisher"), '[[',"count"))

save(res, file = paste0(Outputdir, "/", Outputname,".Rd"))
write.table(resdf, file= paste0(Outputdir,"/", Outputname, ".txt") , row.names=T, col.names=T, sep="\t", quote = F)
warnings()
loginfo("diffdriver finished.")


=======
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

res <- diffdriver(genef, mutf, phenof, drivermapsdir = drivermapsdir,mode=2, k=6,outputdir = outputdir, outputname = outputname)

save(res,file="res.Rd")
>>>>>>> 94530f490d469c3e87bbf5d6ceb53c54b88b60f5


