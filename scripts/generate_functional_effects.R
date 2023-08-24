load("~/SimingLab/jiezhou/diffdriverPlaygournd/EMalgorithm/real/nttype/BRCA/temp/diffDriver_demo_Allergic_disease_TP53.Rd")
load("~/SimingLab/jiezhou/diffdriver/data-raw/annodataFE_OG.Rd")
load("~/SimingLab/jiezhou/diffdriver/data-raw/annodataFE_TSG.Rd")
mutation=function(genenames){
filenames=list.files("~/SimingLab/jiezhou/diffdriverPlaygournd/EMalgorithm/real/nttype/BRCA/temp")
for (i in genenames){
dataname=min(grep(i,filenames))
print(dataname)
}

}
browser()
mutation(names(annodataFE_OG))
