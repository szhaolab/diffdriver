
#' @param BMRmode "regular" or "signature". In regular mode, number of signatures (k) is not used.
#'
#' @return A list
#' @noRd
matrixlistToBMR  <- function(afileinfo, mut, BMRmode = c("signature", "regular"), k=6, outputbase = "."){

  BMRmode <- match.arg(BMRmode)

  # Read annotation files.
  matrixlist <- readmodeldata(afileinfo, yfileinfo = NULL, bmvars, bmmuttype, bmreadinvars, qnvars, functypecodelevel = NULL, qnvarimpute=c(NA))

  # Mutation data to y
  anno <- do.call(rbind, lapply(matrixlist, '[[', 5)) # Extract only chrom, pos, ref, alt and nttype info.
  anno$nttypecode = do.call(rbind, lapply(matrixlist, '[[', 4))
  anno$posidx = do.call(rbind, lapply(matrixlist, function(x) {if (nrow(x[[1]]) >0) data.table::data.table(posidx = 1:nrow(x[[1]]))}))

  id = unique(mut$SampleID)
  nsample = length(id)
    # `ymatrix` is no mut for nttype x sample.
  totalnttype <- length(matrixlist)
  ymatrix = matrix(0, nrow = nsample, ncol = totalnttype)
  rownames(ymatrix) <- id
    # `ysample` is no. mut for each sample.
  ysample <- rep(0, nsample)
  names(ysample) <- id

  for (i in 1:nsample) {
    idata = mut[mut$SampleID==id[i],]
    ni = nrow(idata)
    for (imut in 1:ni) {
      Position = idata$Position[imut]
      Ref = idata$Ref[imut]
      Alt = idata$Alt[imut]
      if (grepl('chr', idata$Chromosome[imut], fixed = TRUE)){
        Chromosome = idata$Chromosome[imut]
      } else{
        Chromosome = paste0("chr",idata$Chromosome[imut])
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

  save(matrixlist, file =  paste0(outputbase,"bmrdebug9-temp.Rd.Rdata"))

  # Infer positional level BMM parameters.
  betabaseline0 <- log(unlist(lapply(lapply(matrixlist,'[[',2), colMeans)))
  min_betabaseline0 <- min(betabaseline0[is.finite(betabaseline0)], na.rm = TRUE) # Replace Inf and NaN with the smallest finite value
  betabaseline0[is.infinite(betabaseline0) | is.nan(betabaseline0)] <- min_betabaseline0

  nbeta <- dim(matrixlist[[1]][[1]])[2] -1
  initpars <- c(betabaseline0, rep(0, nbeta), 0, 2)
  fixstatus <- c(rep(T, totalnttype), rep(F, nbeta), T, F)
  Y_g_s_0 <- data.table::data.table(agg_var = character(), y = numeric(), key = "agg_var")
  Mu_g_s_0 <- data.table::data.table(agg_var = character(), V1 = numeric(), key = "agg_var")
  BMRpars <- optifix(initpars, fixstatus, loglikfn, matrixlist= matrixlist, y_g_s_in=Y_g_s_0, mu_g_s_in=Mu_g_s_0, method = "BFGS", control=list(trace=6, fnscale=-1), hessian=T)
  names(BMRpars$fullpars) <- c(paste("nttype", 1:totalnttype, sep=""), colnames(matrixlist[[1]][[1]])[-1], "beta_f0", "alpha")

  Y_g_s_all <- gene_y(matrixlist)
  Vbeta_s <- BMRpars$fullpars[1: length(BMRpars$fullpars)-1]
  Mu_g_s_all <- gene_mu(Vbeta_s, matrixlist)

  BMRreg <- list("BMRpars" = BMRpars,
                 "Y_g_s_all" = Y_g_s_all,
                 "Mu_g_s_all" = Mu_g_s_all,
                 "nsyn" = sum(ysample),
                 "ysample" = ysample)

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
    colnames(ff)=paste("factor",1:k, sep= "")
    ff=cbind(sigmapping[ ,c("context","alt_allele")], ff)
    sigmtx=ll%*%t(ff[,-c(1,2)]) # rows are samples, column is nttypecode

    # # get positional adjustment for each nt type, based on driverMAPS estimate
    # fixmusdfile <-  system.file("extdata", "colmu_sd_funct78.Rdata", package = "diffdriver")
    # matrixlist <- readmodeldata(afileinfo, yfileinfo = NULL, bmvars, bmmuttype, c("genename", bmvars , "functypecode"), qnvars, functypecodelevel = NULL, qnvarimpute=c(0,0), fixmusd= fixmusdfile) # note normalization of qnvar is approximately right.
    #

    alpha= BMRreg$BMRpars$fullpars[c("alpha")]
    genesp=data.table(genename=BMRreg$Y_g_s_all$agg_var,
                      lambda=(BMRreg$Y_g_s_all$y + alpha)/(BMRreg$Mu_g_s_all$V1 + alpha),
                      key = "genename")
    BMvbeta <-  BMRreg$BMRpars$fullpars[qnvars]
    at <- c()
    for (j in (1:totalnttype)){
      BManno  <- matrixlist[[j]][[1]][, qnvars, with=F]
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

  BMRres <- list("reg" = BMRreg, "sig" = BMRsig)

  BMRfile <- paste0(outputbase,"_BMRres.Rdata")
  save(BMRres, file = BMRfile)

  return(BMRres)
}














