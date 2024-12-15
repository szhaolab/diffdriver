
#'
#' @param annodir The path of the directory that contains annotation files
#' @param mutf The path to the mutation files
#' @param BMRlist The parameters from drivermaps
#' @param target The target positions where the bmr's are computed
#'
#' @return A list
#' @noRd
matrixlistToBMR  <- function(afileinfo, mutf, BMRmode, k=6){
  # Read annotation files.
  matrixlist <- readmodeldata(afileinfo, yfileinfo = NULL, bmvars, bmmuttype, bmreadinvars, qnvars, functypecodelevel = NULL, qnvarimpute=c(-1.8))

  # Read mutation data.
  mutfile <- data.table::fread(mutf, header = T)

  #  Mutation data to y
  chrposmatrixlist <- ddmread(afileinfo, yfileinfo = NULL, selectvars = c("chrom", "start","ref","alt","nttypecode"), bmmuttype, readinvars = c("genename", "chrom", "start", "ref", "alt", "nttypecode", "functypecode"))

  anno <- do.call(rbind, lapply(chrposmatrixlist, '[[', 1)) # Extract only chrom, pos, ref, alt and nttype info.
  anno$nttypecode = do.call(rbind, lapply(chrposmatrixlist, '[[', 4))
  anno$posidx = do.call(lapply(chrposmatrixlist, function(x) 1:nrow(x[[1]])))

  id = unique(mutfile$SampleID)
  nsample = length(id)
    # `ymatrix` is no mut for nttype x sample.
  ymatrix = matrix(0, nrow = nsample, ncol = totalnttype)
  rownames(ymatrix) <- id
    # `ysample` is no. mut for each sample.
  ysample <- rep(0, nsample)
  names(ysample) <- id


  for (i in 1:nsample) {
    idata = mutfile[SampleID==id[i]]
    ni = nrow(idata)
    for (ii in 1:ni) {
      Position = idata$Position[ii]
      Ref = idata$Ref[ii]
      Alt = idata$Alt[ii]
      if (grepl('chr', idata$Chromosome[ii], fixed = TRUE)){
        Chromosome = idata$Chromosome[ii]
      } else{
        Chromosome = paste0("chr",idata$Chromosome[ii])
      }
      mutanno = anno[start==Position & ref==Ref & alt==Alt & chrom==Chromosome, ]
      ntidx = mutanno$nttypecode
      posidx = mutanno$posidx

      if (length(ntidx) == 1){
        ymatrix[i,ntidx]=ymatrix[i,ntidx] + 1
        ysample[i] <-ysample[i] + 1
        matrixlist[[ntidx]][[2]][posidx] <- matrixlist[[ntidx]][[2]][posidx] + 1
      }
    }
  }

  ysample[ysample < 3] <- 3 # if too few syn mutations force it to be 3.
  # TODO: if matrixlist, if one nt got too few variants then add 1.

  # Infer positional level BMM parameters.
  betabaseline0 <-log(apply(rbind(unlist(lapply(lapply(matrixlist,'[[',2), colMeans)),
                                  1/unlist(lapply(lapply(matrixlist,'[[',2), nrow))),2,max))
  nbeta <- dim(matrixlist[[1]][[1]])[2] -1
  initpars <- c(betabaseline0, rep(0, nbeta), 0, 2)
  fixstatus <- c(rep(T, Totalnttype), rep(F, nbeta), T, F)
  Y_g_s_0 <- data.table(agg_var = character(), y = numeric(), key = "agg_var")
  Mu_g_s_0 <- data.table(agg_var = character(), V1 = numeric(), key = "agg_var")
  BMRreg <- optifix(initpars, fixstatus, loglikfn, matrixlist= Matrixlist, y_g_s_in=Y_g_s_0, mu_g_s_in=Mu_g_s_0, method = "BFGS", control=list(trace=6, fnscale=-1), hessian=T)
  names(BMRreg$fullpars) <- c(paste("nttype", 1:Totalnttype, sep=""), colnames(Matrixlist[[1]][[1]])[-1], "beta_f0", "alpha")

  # Adjust for mutational signature in signature mode.
  BMRsig <- NULL
  if (BMRmode == "signature") {

    # Run topic modeling
    # Alternative 1: LDA
    # fit <- topicmodels::LDA(ymatrix, k= 6)
    # posterior(fit)$terms
    # ll <- posterior(fit)$topics
    # ff <- t(posterior(fit)$terms)

    # Alternative 2: fasttopic
    colnames(ymatrix) <- 1: totalnttype
    rownames(ymatrix) <- id
    ymatrix <- ymatrix[rowSums(ymatrix) >= 5, ] # if the sample has too few mutations, use cohort-wise average.

    fit =fastTopics::fit_topic_model(ymatrix, k = k)
    lltemp = fit$L
    fftemp = fit$F
    ff = matrix(min(fftemp)/2, nrow = totalnttype, ncol = k )
    ff[as.numeric(rownames(fftemp)), ] <- fftemp
    ll = matrix(rep(colMeans(lltemp), nsample), nrow = nsample, byrow =T)
    rownames(ll) <- id
    ll[rownames(lltemp),] <- lltemp

    colnames(ll)=paste("weight",1:k, sep = "")
    colnames(ff)=paste("factor",1:k,sep="")
    ff=cbind(sigmapping[,c("context","alt_allele")],ff)
    sigmtx=ll%*%t(ff[,-c(1,2)]) # rows are samples, column is nttypecode

    # get positional adjustment for each nt type, based on driverMAPS estimate
    fixmusdfile <-  system.file("extdata", "colmu_sd_funct78.Rdata", package = "diffdriver")
    matrixlist <- readmodeldata(afileinfo, yfileinfo = NULL, bmvars, bmmuttype, c("genename", bmvars , "functypecode"), qnvars, functypecodelevel = NULL, qnvarimpute=c(0,0), fixmusd= fixmusdfile) # note normalization of qnvar is approximately right.

    BMcol <- c("expr","repl","hic")
    alpha=BMRlist$BMpars$fullpars[c("alpha")]
    genesp=data.table(genename=BMRlist$Y_g_s_all$agg_var,
                      lambda=(BMRlist$Y_g_s_all$y + alpha)/(BMRlist$Mu_g_s_all$V1 + alpha),
                      key = "genename")
    BMvbeta <- BMRlist$BMpars$fullpars[c("expr","repl","hic")]
    at <- c()
    for (j in (1:totalnttype)){
      BManno  <- matrixlist[[j]][[1]][, BMcol, with=F]
      genename <- matrixlist[[j]][[3]]
      at[j] <- sum(exp(as.matrix(BManno) %*% BMvbeta) * as.matrix(genesp[genename,"lambda"]), na.rm = T)
    }

    # some nttype doesn't have silent mutations.
    at[at==0] <- max(at)

    BMRsig <- list(sigmtx = sigmtx,
                   at = at,
                   BMvbeta = BMvbeta,
                   genesp = genesp,
                   ll = ll,
                   ff = ff,
                  ysample =ysample,
                   fit = fit)
  }

  bmrres <- list("reg" = BMRreg, "sig" = BMRsig)

  return(bmrres)
}














