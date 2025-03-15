
prep_positional_data <- function(gene, afileinfo, BMRreg = NULL, output_prefix = NULL, output_dir = "."){
  gene <- data.frame(gene)
  fixmusdfile <- system.file("extdata", "colmu_sd_funct78.Rdata", package = "diffdriver")
  genef <- tempfile(output_prefix, tmpdir = output_dir )
  write.table(gene, file = genef, col.names = F, quote = F, row.names = F)
  matrixlist <- readmodeldata(afileinfo, yfileinfo = NULL, c(bmvars,funcvars), funcvmuttype, readinvars , qnvars, functypecodelevel, qnvarimpute=c(-1.8), cvarimpute = 0, genesubset=genef, fixmusd= fixmusdfile)
  chrposmatrixlist <- ddmread(afileinfo, yfileinfo = NULL, c("chrom", "start","ref","alt","nttypecode"), funcvmuttype, c("genename", "chrom", "start", "ref", "alt", "functypecode", "ssp", "nttypecode"), genesubset=genef)
  file.remove(genef)

  # ------ debug------------------
  # save(matrixlist, chrposmatrixlist, file = '/dartfs-hpc/rc/home/m/f0052zm/temp/output/debug_glmdt_96_sig.Rd')
  # load('/dartfs-hpc/rc/home/m/f0052zm/temp/output/debug_glmdt_96_sig.Rd')
  totalnttype <- length(matrixlist)
  for (t in 1:  totalnttype){
    b1 <- which(sapply(matrixlist[[t]][[3]], grepl, pattern="[;,|]"))
    matrixlist[[t]][[3]][b1,] <- unlist(lapply(sapply(matrixlist[[t]][[3]][b1,], strsplit, split="[;,|]"), function(x) intersect(x,gene[,1])[1]))
    b2 <- which(sapply(chrposmatrixlist[[t]][[3]], grepl, pattern="[;,|]"))
    chrposmatrixlist[[t]][[3]][b2,] <- unlist(lapply(sapply(matrixlist[[t]][[3]][b2,], strsplit, split="[;,|]"), function(x) intersect(x,gene[,1])[1]))
  }

  matrixlisttemp <- copy(matrixlist)
  chrposmatrixlisttemp <- copy(chrposmatrixlist)
  y_g_s <- BMRreg$Y_g_s_all[gene][,1:2, with=F]
  y_g_s[is.na(y_g_s)] <- 0
  mu_g_s <- BMRreg$Mu_g_s_all[gene][,1:2, with=F]
  mu_g_s[is.na(mu_g_s)] <- 0

  glmdtall <- matrixlistToGLM(matrixlisttemp, chrposmatrixlisttemp, BMRreg$BMRpars, mu_g_s, y_g_s, fixpars=NULL)
  rm(matrixlisttemp,chrposmatrixlisttemp,matrixlist,chrposmatrixlist); gc()

  # ------ debug -------------------
  # save(glmdtall, file = '/dartfs-hpc/rc/home/m/f0052zm/temp/output/debug_glmdt_96_sig.Rd')

  # functional annotation (fanno) :data.table, each row has functional annotation for each possible mutation
  fanno <- glmdtall[[1]]

  # row index (ri): chr pos ref alt
  ri <- glmdtall[[2]][,.(chrom,genename,start,ref,alt,nttypecode)]

  return(list("fanno" = fanno,
              "ri"  = ri,
              "BMRbaseline" = glmdtall[[1]]$baseline))

}


# mutations (mut): data.table, with columns Chromosome, Position, Ref, Alt, SampleID
prep_pheno_mut_data <- function(mut, pheno){
  mut <- data.table::data.table(mut)
  if (!grepl('chr', mut$Chromosome[1], fixed = T)) {mut$Chromosome <- paste0("chr",mut$Chromosome)}

  # sample annotation (pheno):data.table, with columns BMR label, No. syn and phenotype.
  pheno <- data.table::data.table(pheno)

  shared=intersect(mut$SampleID,pheno$SampleID)
  if (length(shared) < 10){
    stop("too few samples with both phenotype and mutation data.")
  }
  index1=which(pheno$SampleID %in% shared)
  index2=which(mut$SampleID %in% shared)

  pheno = pheno[index1,]
  mut<- mut[index2,]

  # column index (ci): sampleID
  ci <- pheno[,"SampleID"]
  ci[,"cidx" := 1:dim(pheno)[1]]

  return(list("pheno" = pheno,
              "mut"   = mut,
              "ci"    = ci))
}





