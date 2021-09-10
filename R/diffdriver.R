# TODO: mutation hotspot, quantitative phenotype
# Notes: syn data prepared from direct join of mutation list to annodata is different from used in driverMAPS BMR estimation: 1. 5% syn are ssp are included. 2. BMR estimation from old driverMAPS run has different annodata. 3. genes without expr/hic/rep were removed in BMR driverMAPS estimation.
#' @title Run diffDriver given input files
#' @description This is the function to run diffDriver. We first need to set up: run driverMAPS for groups with potential different BMR, assign BMR labels for each sample. then BMR for each sample will be scaled based on driverMAPS results.
#' @param genef file for name of genes to be included in the analysis
#' @param mutf mutation list file, use the driverMAPS mutation input format
#' @param phenof phenptype file, SampleID <tab> Phenotype <tab> Nsyn. nsyn is number of syn mutations in this sample.
#' @import Matrix data.table
#' @export
diffDriver <- function(genef, mutf, phenof){
  # prep data
  # column index (ci): sampleID
  # row index (ri): chr pos ref alt
  #
  # 1. background mutation rate (bmrdt): data.table, each column correponds to one BMR label
  # 2. functional annotation (fanno):data.table, each row has functional annotation for each possible mutation
  # 3. sample annotation (canno):data.table, with columns BMR label, No. syn and phenotype.

  allg <- read.table(genef, stringsAsFactors = F)[,1]
  matrixlist <- readmodeldata(Afileinfo, Yfileinfo, c(BMvars,Funcvars), Funcvmuttype, Readinvars , Qnvars, Functypecodelevel,qnvarimpute=c(0,0), cvarimpute = 0, genesubset=genef, fixmusd=Fixmusdfile)
  chrposmatrixlist <- ddmread(Afileinfo, Yfileinfo, c("chrom", "start","ref","alt"), Funcvmuttype, c("genename", "chrom", "start", "ref", "alt", "functypecode", "ssp", "nttypecode"), genesubset=genef)

  for (t in 1:length(matrixlist)){
    b1 <- which(sapply(matrixlist[[t]][[3]], grepl, pattern="[;,|]"))
    matrixlist[[t]][[3]][b1,] <- unlist(lapply(sapply(matrixlist[[t]][[3]][b1,], strsplit, split="[;,|]"), function(x) intersect(x,allg)[1]))
    b2 <- which(sapply(chrposmatrixlist[[t]][[3]], grepl, pattern="[;,|]"))
    chrposmatrixlist[[t]][[3]][b2,] <- unlist(lapply(sapply(matrixlist[[t]][[3]][b2,], strsplit, split="[;,|]"), function(x) intersect(x,allg)[1]))
  }

  bmrdt <- data.table()

  for (label in names(BMRlist)) {
    matrixlisttemp <- copy(matrixlist)
    chrposmatrixlisttemp <- copy(chrposmatrixlist)
    y_g_s <- BMRlist[[label]]$Y_g_s_all[allg][,1:2, with=F]
    y_g_s[is.na(y_g_s)] <- 0
    mu_g_s <- BMRlist[[label]]$Mu_g_s_all[allg][,1:2, with=F]
    mu_g_s[is.na(mu_g_s)] <- 0
    glmdtall <- matrixlistToGLM(matrixlisttemp, chrposmatrixlisttemp, BMRlist[[label]]$BMpars, mu_g_s, y_g_s, fixpars=NULL)
    bmrdt[,eval(label):= glmdtall[[1]]$baseline]
    rm(matrixlisttemp,chrposmatrixlisttemp); gc()
  }

  fanno <- glmdtall[[1]]

  ri <- glmdtall[[2]]

  canno <- fread(phenof, header = T)
  ci <- canno[,"SampleID"]
  ci[,"cidx" := 1:dim(canno)[1]]
  for (l in names(BMRlist)) {
    BMRlist[[l]][["nsyn"]] <- sum(canno$Nsyn[canno$BMRlabel == l])
  }

  muts <- fread(mutf, header = T)
  if (!grepl('chr', muts$Chromosome[1], fixed = T)) {muts$Chromosome <- paste0("chr",muts$Chromosome)}

  # split based on gene
  bmrallg <- split(bmrdt, ri$genename)
  riallg <- split(ri,ri$genename)
  fannoallg <- split(fanno,ri$genename)

  # run diffdriver for each gene
  res <- list()
  for (g in names(bmrallg)) {
    print(paste0("Start to process gene: ", g))
    rig <- riallg[[g]]
    rig$ridx <- 1:dim(rig)[1]
    muti <- na.omit(ci[rig[muts, on = c("chrom"= "Chromosome", "start" = "Position",  "ref" = "Ref",  "alt"= "Alt")], on = "SampleID"])
    mutmtx <- sparseMatrix(i = muti$ridx, j = muti$cidx, dims = c(max(rig$ridx), max(ci$cidx)))
    if (sum(mutmtx) ==0) {next}

    bmrsc <- log(canno$Nsyn/unlist(lapply(BMRlist[canno$BMRlabel],'[[',"nsyn"))) # normalize BMR for each sample
    # bmrmtx is matrix, rows are positions, columns are for each samples
    bmrmtx <-sweep(as.matrix(bmrallg[[g]][,canno$BMRlabel,with=F]), 2, bmrsc, "+")

    betaf <- Fpars[[g]][names(Fpars[[g]]) != "beta_f0"]
    betaf0 <- Fpars[[g]]["beta_f0"]
    if (is.null(betaf)){
      betaf <-  Fpars[["TP53"]][names(Fpars[["TP53"]]) != "beta_f0"]
      betaf0 <-  Fpars[["TP53"]]["beta_f0"]
    } # if OG/TSG unknown, use TSG parameters.
    ganno <- fannoallg[[g]]
    fe <- as.matrix(ganno[ ,names(betaf), with =F]) %*% betaf + betaf0

    resg <- list()
    resg[["dd"]] <- ddmodel(mutmtx, canno$Phenotype, bmrmtx, fe[,1])
    resg[["mlr"]] <- mlr(mutmtx, canno$Phenotype)
    resg[["mlr.v2"]] <- mlr.v2(mutmtx, canno$Phenotype, canno$Nsyn)
    resg[["fisher"]] <- genefisher(mutmtx, canno$Phenotype)
    resg[["binom"]] <- genebinom(mutmtx, canno$Phenotype)
    resg[["lr"]] <- genelr(mutmtx, canno$Phenotype)
    res[[g]] <- resg

    save(mutmtx, canno, bmrmtx, fe, ganno, betaf, betaf0, resg, file=paste0(paste0(Outputbase,".", g, ".Rd")))
    setEPS()
    postscript(file=paste0(Outputbase,".", g, "mut_status.eps"), width=9, height=4)
    plot_mut(mutmtx, canno, bmrmtx, ganno)
    dev.off()
  }

  return(res)
}