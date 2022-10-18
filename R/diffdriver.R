# TODO: mutation hotspot, quantitative phenotype
# Notes: syn data prepared from direct join of mutation list to annodata is different from used in driverMAPS BMR estimation: 1. 5% syn are ssp are included. 2. BMR estimation from old driverMAPS run has different annodata. 3. genes without expr/hic/rep were removed in BMR driverMAPS estimation.
#' @title Detect the relationship between a specified phenotype and mutation
#' @description This is the function to run diffDriver. We first need to set up: run driverMAPS for groups with potential different BMR, assign BMR labels for each sample. then BMR for each sample will be scaled based on driverMAPS results.
#' @param genef Link for files of name of genes to be included in the analysis
#' @param mutf Link for mutation file, use the driverMAPS mutation input format
#' @param phenof Link for phenoptype file, SampleID <tab> Phenotype <tab> Nsyn. nsyn is number of syn mutations in this sample.
#' @param mode Mode indicating the type of nucleotide change. The case of mode=1 corresponds to
#' regular change while mode=2 the mutation signature case.
#' @param j The index of phenotype
#' @param k The specified number of components in NMF
#' @param model The type of annotation. The case of mode=1 is for regular annotation while the case of mode=2 is for signature annotation file
#' @return A list
#' @export
diffdriver <- function(genef, mutf, phenof,j, hotf=NULL, drivermapsdir,k=5,mode=1, outputdir =".", outputname = "diffdriver_results"){
if (mode==1){
  res <- diffdriver_reg(genef, mutf, phenof, j, hotf, drivermapsdir = drivermapsdir, outputdir = outputdir, outputname = outputname)
}else{
  res <- diffdriver_sig(genef, mutf, phenof, j, hotf, drivermapsdir = drivermapsdir,k=k, outputdir = outputdir, outputname = outputname)
}
  return(res)
}
