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
#'   \item{Phenotype1}{<dbl>}
#'   \item{Phenotype2}{<dbl>}
#'   \item{...}{...}
#' }
#' Example:
#' \preformatted{
#' SampleID SmokingCessation BMI
#' TCGA-N5-A4R8-01A-11D-A28R-08 0.5319630 20.0
#' TCGA-N5-A4RD-01A-11D-A28R-08 0.0448991 24.4
#' }
#'
#' @param anno_dir The path to the directory with all the annotation files.
#'  Please download from Zenodo.The default is current folder
#'
#' @param totalnttype either 9 or 96. Will look for annotation files
#' anno9_ntypexxx_annodata.txt when totalnttype is 9 or anno96_ntypexxx_annodata
#' when totalnttype is 96.
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
                     totalnttype = 96,
                     BMRmode = c("signature", "regular"),
                     output_dir =".",
                     output_prefix = "diffdriver_results"){

  # ------- SET UP ----------
  BMRmode <- match.arg(BMRmode)

  afileinfo <- list(file = file.path(anno_dir, paste0("anno", totalnttype, "_nttypeXXX_annodata.txt")),
                      header = aheader,
                      coltype = acoltype,
                      totalntype = totalnttype)

  dir.create(output_dir)

  outputbase <- paste0(output_dir, "/", output_prefix)

  # Read in silent mutation data and infer BMR
  # note silent mutations from all samples in `mut` are used in getting
  # BMR parameter estimates, not just the ones with phenotype info.
  print("Infer parameters in background mutation rate model")

  # BMRres <- matrixlistToBMR(afileinfo, mut = mut, BMRmode, k=k, outputbase)

  # ------ debug------------------
  load("~/temp/output/debug_sig_BMRres.Rdata")

  BMRreg <- BMRres[[1]]

  if (BMRmode == "signature"){
    bmrsig <- BMRres[[2]]
  }

  # ------- Prep data for target genes----------
  print("Start to prepare input data for target genes ...")

  rdata <- prep_positional_data(gene, afileinfo, BMRreg = BMRreg, output_prefix = output_prefix, output_dir =  output_dir)
  fanno <-  rdata$fanno
  ri <- rdata$ri

  nsyndt <- data.table::data.table("SampleID" = names(BMRreg[["ysample"]]), "Nsyn" = BMRreg[["ysample"]])
  cdata <- prep_pheno_mut_data(mut, pheno, add_dt = nsyndt)
  mut <-cdata$mut
  ci <- cdata$ci
  pheno <- cdata$pheno

  # add hotspots
  hotspots <- mut2hotspot(mut)

  ## ------- Get BMR (log scale mu_ij) for target genes -----------------------------

  if (BMRmode == "regular"){
    bmrdt <- data.table()
    bmrdt[,"BMR":= rdata$BMRbaseline]
    bmrsc <- log(pheno$Nsyn/BMRreg$nsyn)
    bmrmtx_uni <- as.matrix(bmrdt[,rep(1,nrow(pheno)), with = F])
    bmrmtx <-as.data.table(sweep(bmrmtx_uni, 2, bmrsc, "+"))
    bmrallg <- split(bmrmtx, ri$genename)
  } else{
    bmr=data.table()
    for (i in 1:nrow(ci)){
      id <- ci[i,"SampleID"][[1]]
      glambda <- bmrsig$genesp[ri$genename, lambda]
      glambda[which(is.na(glambda))] <- 1
      mui <- bmrsig$ysample[id] * bmrsig$sigmtx[id,ri$nttypecode] *
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
    muti <- na.omit(ci[rig[mut, on = c("chrom"= "Chromosome", "start" = "Position",  "ref" = "Ref",  "alt"= "Alt")], on = "SampleID"])
    mutmtx <- Matrix::sparseMatrix(i = muti$ridx, j = muti$cidx, dims = c(max(rig$ridx), max(ci$cidx)))

    hotg= na.omit(rig[hotspots,on=c("chrom"="Chromosome","start" = "Position")])
    hotmat=rep(0,nrow(rig))
    hotmat[hotg$ridx]=1
    if (sum(mutmtx) == 0) {
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
    e=pheno[[cdata$phename]]
    resg[["dd"]] <- ddmodel(mutmtx, e, bmrmtx, fe[,1], label=label)
    ## resg[["dd_nl"]] <- ddmodel_nl(mutmtx, e, bmrmtx, fe[,1])
    resg[["mlr"]] <- mlr(mutmtx, e)
    resg[["mlr.v2"]] <- mlr.v2(mutmtx, e, pheno$Nsyn)
    e_binary=ifelse(e>mean(e),1,0)
    resg[["fisher"]] <- genefisher(mutmtx, e_binary)
    resg[["binom"]] <- genebinom(mutmtx, e_binary)
    resg[["lr"]] <- genelr(mutmtx, e_binary)
    res[[g]] <- resg
  }


  if (BMRmode == "signature") {
    fastopicfit <- bmrsig$fit
    save(fastopicfit, e, bmrallg, fannoallg, ci, riallg, res, file=paste0(outputbase, "_" , cdata$phename, "_resdd.Rd"))
  } else {
    save(e, bmrallg, fannoallg, ci, riallg, res, file=paste0(outputbase, "_" , cdata$phename, "_resdd.Rd"))
  }

  meth <- c("dd", "mlr", "mlr.v2", "fisher", "binom", "lr")
  resdf <- data.frame(lapply(meth, function(x)unlist(lapply(lapply(res,'[[',x),'[[','pvalue'))))
  colnames(resdf) <- paste0(meth,".p")
  resdf[ , paste0(meth,".fdr")] <- apply(resdf,2,p.adjust, method = "fdr")
  resdf[,c("mut.E1", "mut.E0", "E1", "E0")] <- do.call(rbind, lapply(lapply(res, '[[', "fisher"), '[[',"count"))
  write.table(resdf, file = paste0(outputbase,"_", cdata$phename, "_resdd.txt"), quote = F, col.names = T, row.names = T)

  print("Finished.")

  return(resdf)
}
