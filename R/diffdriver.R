# TODO: mutation hotspot, quantitative phenotype
# Notes: syn data prepared from direct join of mutation list to annodata is different from used in driverMAPS BMR estimation: 1. 5% syn are ssp are included. 2. BMR estimation from old driverMAPS run has different annodata. 3. genes without expr/hic/rep were removed in BMR driverMAPS estimation.
#' @title Run diffDriver given input files
#' @description This is the function to run diffDriver. We first need to set up: run driverMAPS for groups with potential different BMR, assign BMR labels for each sample. then BMR for each sample will be scaled based on driverMAPS results.
#' @param genef file for name of genes to be included in the analysis
#' @param mutf mutation list file, use the driverMAPS mutation input format
#' @param phenof phenptype file, SampleID <tab> Phenotype <tab> Nsyn. nsyn is number of syn mutations in this sample.
#' @param mode Mode indicating the type of nucleotide change. The case of mode=1 corresponds to
#' regular change while mode=2 the mutation signature case.
#' @import Matrix data.table
#' @export
diffdriver <- function(genef, mutf, phenof, drivermapsdir,mode=1, outputdir =".", outputname = "diffdriver_results"){
if (mode==1){
  res <- diffdriver_reg(genef, mutf, phenof, drivermapsdir = drivermapsdir, outputdir = outputdir, outputname = outputname)
}else{
  res <- diffdriver_sig(genef, mutf, phenof, drivermapsdir = drivermapsdir, outputdir = outputdir, outputname = outputname)
}
  return(res)
}
