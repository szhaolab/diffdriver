library("diffdriver")
library("brglm")
library("fastTopics")
library("data.table")
drivermapsdir <- "~/SimingLab/library/diffdriver_anno/"
outputname <- "temp/diffDriver_demo"
i1=17
i2=9
typee=list.files("/dartfs/rc/lab/S/Szhao/diffDriver/data_run_prs/tumor_specific_input/")
type=typee[i1]
outputdir <- paste0("~/SimingLab/jiezhou/diffdriverPlaygournd/EMalgorithm/real/signature/",type)
mutationdir="/dartfs/rc/lab/S/Szhao/diffDriver/data_run_prs/tumor_specific_input/"
genef = paste0(mutationdir,type,"/",type,"genes.txt")
mutf = paste0(mutationdir,type,"/",type,"_mutations.txt")
phenof =paste0(mutationdir,type,"/",type,"_PRS_46phenotype.txt")
pheName=colnames(fread(phenof, header = "auto"))
hotf=paste0(mutationdir,type,"/",type,"_hotpositions.txt")
res <- diffdriver(genef, mutf, phenof,j=i2+4,hotf, drivermapsdir = drivermapsdir,k=6,mode=2, outputdir = outputdir, outputname = outputname)
print(paste0(type,":",i2+4))
save(res,file=paste0(outputdir,"/Results_",type,"_",pheName[i2+4],".Rd"))
q("no")
