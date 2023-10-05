library("diffdriver")
library("data.table")
library("brglm")
drivermapsdir <- "/dartfs/rc/lab/S/Szhao/jiezhou/driverMAPS/data/"
outputname <- "temp/diffDriver_demo"
i1=17
i2=9
typee=list.files("/dartfs/rc/lab/S/Szhao/diffDriver/data_run_prs/tumor_specific_input/")
type=typee[i1]
outputdir <- paste0("~/SimingLab/jiezhou/diffdriverPlaygournd/EMalgorithm/real/nttype/",type)
mutationdir="/dartfs/rc/lab/S/Szhao/diffDriver/data_run_prs/tumor_specific_input/"
genef = paste0(mutationdir,type,"/",type,"genes.txt")
mutf = paste0(mutationdir,type,"/",type,"_mutations.txt")
phenof =paste0(mutationdir,type,"/",type,"_PRS_46phenotype.txt")
bmrf=paste0(mutationdir,type,"/","BMRlist.Rd")

pheName=colnames(data.table::fread(phenof, header = "auto"))
hotf=paste0(mutationdir,type,"/",type,"_hotpositions.txt")
res <- diffdriver(genef, mutf, phenof,bmrf=bmrf,j=i2+4,hotf, drivermapsdir = drivermapsdir,k=6,
                  mode=1, outputdir = outputdir, outputname = outputname)
save(res,file=paste0(outputdir,"/Results_",type,"_",pheName[i2+4], ".Rd"))
q("no")
