#' @title Run diffDriver with Input Files
#' @description This function runs diffDriver.
#' @param gene A vector of genes to be included in the analysis.
#' @param mut A data frame containing all somatic mutations from the cohort. The format is:
#' \describe{
#'   \item{Chromosome}{<int>}
#'   \item{Position}{<int>}
#'   \item{Ref}{<chr>}
#'   \item{Alt}{<chr>}
#'   \item{SampleID}{<chr>}
#' }
#' Example:
#' \preformatted{
#' Chromosome Position Ref Alt SampleID
#' 1 19 55653236 C T TCGA-N6-A4VE-01A-11D-A28R-08
#' }
#' @param pheno A data frame containing sample phenotypes. The format is:
#' #' \describe{
#'   \item{SampleID}{<chr>}
#'   \item{Phenotype}{<dbl>}
#' }
#' Example:
#' \preformatted{
#' SampleID SmokingCessation
#' 1 TCGA-N5-A4R8-01A-11D-A28R-08 0.5319630
#' 2 TCGA-N5-A4RD-01A-11D-A28R-08 0.0448991
#' }
#'
#' @param anno_dir The path to the directory with all the annotation files.
#'  Please download from Zenodo.The default is current folder
#'
#' @param k The number of topics used in modeling background mutation rate.
#'  The default is 6.
#'
#' @param BMRmode There are two modes to run diffdriver. One is "signature",
#' this will model individual level BMR, this is the default. The second one is
#' "regular", this assumes BMR is the same across individuals, only models
#' position-level difference.
#'
#' @param output_dir The path to output directory
#'
#' @param output_prefix The prefix being added to the output file names.
#'
#' @import Matrix
#' @import data.table
#' @export

diffdriver= function(gene,
                     mut,
                     pheno,
                     anno_dir = ".",
                     k = 6,
                     BMRmode = c("signature", "regular"),
                     output_dir =".",
                     output_prefix = "diffdriver_results"){

  # ------- SET UP ----------
  BMRmode <- match.arg(BMRmode)
  if (!(BMRmode %in% c("regular", "signature"))){
    stop("Unknown BMR mode, options: signature, regular")
  }

  afileinfo <- list(file = file.path(annodir, "anno96_nttypeXXX_annodata.txt"),
                      header = aheader,
                      coltype = acoltype)
  totalnttype <<- 96

  dir.create(outputdir)

  outputbase <<- paste0(outputdir, "/", outputname)

  # Read in silent mutation data and infer BMR
  # note silent mutations from all samples in `mut` are used in getting
  # BMR parameter estimates, not just the ones with phenotype info.
  print("Infer parameters in background mutation rate model, adjusting for \
        inter individual mutational signature difference ...")
  BMRres <- matrixlistToBMR(afileinfo, mut, BMRmode, k=k)

  BMRlist <- BMRres[[1]]

  if (BMRmode == "signature"){
    bmrsig <- BMRres[[2]]
  }


  # ------- Read in data for target genes----------
  print("Start to read in data for target genes ...")
  allg <- read.table(genef, stringsAsFactors = F)[,1]


  fixmusdfile <- system.file("extdata", "colmu_sd_funct78.Rdata", package = "diffdriver")
  matrixlist <- readmodeldata(afileinfo, yfileinfo = NULL, c(bmvars,funcvars), funcvmuttype, readinvars , qnvars, functypecodelevel, qnvarimpute=c(-1.8), cvarimpute = 0, genesubset=genef, fixmusd= fixmusdfile)
  chrposmatrixlist <- ddmread(afileinfo, yfileinfo = NULL, c("chrom", "start","ref","alt","nttypecode"), funcvmuttype, c("genename", "chrom", "start", "ref", "alt", "functypecode", "ssp", "nttypecode"), genesubset=genef)


  for (t in 1:length(matrixlist)){
    b1 <- which(sapply(matrixlist[[t]][[3]], grepl, pattern="[;,|]"))
    matrixlist[[t]][[3]][b1,] <- unlist(lapply(sapply(matrixlist[[t]][[3]][b1,], strsplit, split="[;,|]"), function(x) intersect(x,allg)[1]))
    b2 <- which(sapply(chrposmatrixlist[[t]][[3]], grepl, pattern="[;,|]"))
    chrposmatrixlist[[t]][[3]][b2,] <- unlist(lapply(sapply(matrixlist[[t]][[3]][b2,], strsplit, split="[;,|]"), function(x) intersect(x,allg)[1]))
  }

  matrixlisttemp <- copy(matrixlist)
  chrposmatrixlisttemp <- copy(chrposmatrixlist)
  y_g_s <- BMRlist$Y_g_s_all[allg][,1:2, with=F]
  y_g_s[is.na(y_g_s)] <- 0
  mu_g_s <- BMRlist$Mu_g_s_all[allg][,1:2, with=F]
  mu_g_s[is.na(mu_g_s)] <- 0
  glmdtall <- matrixlistToGLM(matrixlisttemp, chrposmatrixlisttemp, BMRlist$BMpars, mu_g_s, y_g_s, fixpars=NULL)
  rm(matrixlisttemp,chrposmatrixlisttemp,matrixlist,chrposmatrixlist); gc()

  # functional annotation (fanno) :data.table, each row has functional annotation for each possible mutation
  fanno <- glmdtall[[1]]

  # row index (ri): chr pos ref alt
  ri <- glmdtall[[2]][,.(chrom,genename,start,ref,alt,nttypecode)]

  # mutations (muts): data.table, with columns Chromosome, Position, Ref, Alt, SampleID
  muts0 <- data.table::fread(mutf, header = T)
  if (!grepl('chr', muts0$Chromosome[1], fixed = T)) {muts0$Chromosome <- paste0("chr",muts0$Chromosome)}

  # sample annotation (canno):data.table, with columns BMR label, No. syn and phenotype.
  canno0 <- data.table::fread(phenof, header = "auto")
  shared=intersect(muts0$SampleID,canno0$SampleID)
  if (length(shared) < 10){
    stop("two few samples with both phenotype and mutation data.")
  }
  index1=which(canno0$SampleID %in% shared)
  index2=which(muts0$SampleID %in% shared)

  canno = canno0[index1,]
  muts<- muts0[index2,]

  # column index (ci): sampleID
  ci <- canno[,"SampleID"]
  ci[,"cidx" := 1:dim(canno)[1]]

  # hotspot: TODO add a function mutf2hotspot
  # hotspots <- mutf2hotspot(mutf)
  hotspots=read.table(file = hotf, header =T)

  ## ------- Get BMR (log scale mu_ij) for target genes -----------------------------

  if (BMRmode == "regular"){
    bmrdt <- data.table()
    bmrdt[,"BMR":= glmdtall[[1]]$baseline]
    BMRlist[["nsyn"]] <- sum(canno$Nsyn)
    bmrsc <- log(canno$Nsyn/BMRlist$nsyn)
    bmrmtx_uni <- as.matrix(bmrdt[,rep(1,nrow(canno)), with = F])
    bmrmtx <-as.data.table(sweep(bmrmtx_uni, 2, bmrsc, "+"))
    bmrallg <- split(bmrmtx, ri$genename)
  } else{
    bmr=data.table()
    for (i in 1:nrow(ci)){
      id <- ci[i,"SampleID"][[1]]
      glambda <- bmrsig$genesp[ri$genename, lambda]
      glambda[which(is.na(glambda))] <- 1
      mui <- bmrsig$yn[id] * bmrsig$sigmtx[id,ri$nttypecode] *
        exp(as.matrix(fanno[, c("expr","repl","hic")]) %*% bmrsig$BMvbeta) *
        glambda/bmrsig$at[ri$nttypecode]
      bmr=cbind(bmr,log(mui))
    }
    if (any(is.na(bmr))) {stop("bmr missing")}
    bmrallg <- split(bmr,ri$genename)
  }

  riallg <- split(ri,ri$genename)
  fannoallg <- split(fanno,ri$genename)

  rm(ri,fanno)

  # run diffdriver for each gene
  res <- list()
  for (g in names(bmrallg)) {
    print(paste0("Start to process gene: ", g))
    rig <- riallg[[g]]
    rig$ridx <- 1:dim(rig)[1]
    muti <- na.omit(ci[rig[muts, on = c("chrom"= "Chromosome", "start" = "Position",  "ref" = "Ref",  "alt"= "Alt")], on = "SampleID"])
    mutmtx <- Matrix::sparseMatrix(i = muti$ridx, j = muti$cidx, dims = c(max(rig$ridx), max(ci$cidx)))

    hotg= na.omit(rig[hotspots,on=c("chrom"="chrom","start"="start")])
    hotmat=rep(0,nrow(rig))
    hotmat[hotg$ridx]=1
    if (sum(mutmtx) ==0) {
      next
    }

    bmrmtx= bmrallg[[g]]
    ganno <- fannoallg[[g]]

    if (g %in% OGs[,1]){
      betaf <- OGpars[names(OGpars) != "beta_f0"]
      betaf0 <- OGpars["beta_f0"]
      fe <-  as.matrix(ganno[ ,names(betaf), with =F]) %*% betaf + hotmat*hmm[9] + betaf0
    } else {
      betaf <- TSGpars[names(TSGpars) != "beta_f0"]
      betaf0 <- TSGpars["beta_f0"]
      fe <-  as.matrix(ganno[ ,names(betaf), with =F]) %*% betaf + betaf0
    } # if OG/TSG unknown, use TSG parameters.

    label=factor(1:nrow(bmrmtx))

    resg <- list()
    e=canno[[j]]
    phename= gsub(" ", "_", colnames(canno)[j])
    resg[["dd"]] <- ddmodel(mutmtx, e, bmrmtx, fe[,1], label=label)
    ## resg[["dd_nl"]] <- ddmodel_nl(mutmtx, e, bmrmtx, fe[,1])
    resg[["mlr"]] <- mlr(mutmtx, e)
    resg[["mlr.v2"]] <- mlr.v2(mutmtx, e, canno$Nsyn)
    e_binary=ifelse(e>mean(e),1,0)
    resg[["fisher"]] <- genefisher(mutmtx, e_binary)
    resg[["binom"]] <- genebinom(mutmtx, e_binary)
    resg[["lr"]] <- genelr(mutmtx, e_binary)
    res[[g]] <- resg
  }


  if (BMRmode == "signature") {
    fastopicfit <- bmrsig$fit
    save( fastopicfit, e, bmrallg, fannoallg, ci, riallg, res, file=paste0(outputbase, "_" , phename, "_resdd.Rd"))
  } else {
    save(e, bmrallg, fannoallg, ci, riallg, res, file=paste0(outputbase, "_" , phename, "_resdd.Rd"))
  }

  meth <- c("dd", "mlr", "mlr.v2", "fisher", "binom", "lr")
  resdf <- data.frame(lapply(meth, function(x)unlist(lapply(lapply(res,'[[',x),'[[','pvalue'))))
  colnames(resdf) <- paste0(meth,".p")
  resdf[ , paste0(meth,".fdr")] <- apply(resdf,2,p.adjust, method = "fdr")
  resdf[,c("mut.E1", "mut.E0", "E1", "E0")] <- do.call(rbind, lapply(lapply(res, '[[', "fisher"), '[[',"count"))
  write.table(resdf, file = paste0(outputbase,"_",phename, "_resdd.txt"))

  print("Finished.")

  return(resdf)
}
