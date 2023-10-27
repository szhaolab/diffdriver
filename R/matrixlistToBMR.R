
#'
#' @param annodir The path of the directory that contains annotation files
#' @param mutf The path to the mutation files
#' @param BMRlist The parameters from drivermaps
#' @param target The target positions where the bmr's are computed
#'
#' @return A list
#' @noRd
matrixlistToBMR  <- function(afileinfo, mutf, BMRlist, k=6){
  # Read annotation files.
  chrposmatrixlist <- ddmread(afileinfo, yfileinfo = NULL, selectvars = c("chrom", "start","ref","alt","nttypecode"), bmmuttype, readinvars = c("genename", "chrom", "start", "ref", "alt", "nttypecode", "functypecode"))

  anno <- do.call(rbind, lapply(chrposmatrixlist, '[[', 1))
  anno$nttypecode = do.call(rbind, lapply(chrposmatrixlist, '[[', 4))

  # Read mutation data
  mutfile <- data.table::fread(mutf, header = T)
  id = unique(mutfile$SampleID)
  nsample = length(id)
  ymatrix = matrix(0, nrow = nsample, ncol = totalnttype)
  rownames(ymatrix) <- id
  yn <- rep(0, nsample) # nonsyn for each sample
  names(yn) <- id
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
      index=anno[start==Position & ref==Ref & alt==Alt & chrom==Chromosome ,]$nttypecode
      if (length(index) == 1){
        ymatrix[i,index]=ymatrix[i,index] + 1
        yn[i] <- yn[i] + 1
      }
    }
  }

  yn[yn < 3] <- 3 # if too few syn mutations force it to be 3.

  # Run topic modeling
  # Alternative: LDA
  # fit <- topicmodels::LDA(ymatrix, k= 6)
  # posterior(fit)$terms
  # ll <- posterior(fit)$topics
  # ff <- t(posterior(fit)$terms)

  # Alternative: fasttopic
  fit =fastTopics::fit_poisson_nmf(ymatrix,k = k,numiter = 100)
  ll = fit$L
  ff = fit$F

  colnames(ll)=paste("weight",1:k, sep = "")
  colnames(ff)=paste("factor",1:k,sep="")
  ff=cbind(sigmapping[,c("context","alt_allele")],ff)
  sigmtx=ll%*%t(ff[,-c(1,2)])

  rownames(sigmtx) <- id # rows are samples, column is nttypecode

  # get positional adjustment for each nt type, based on driverMAPS estimate
  fixmusdfile <-  system.file("extdata", "colmu_sd_funct78.Rdata", package = "diffdriver")
  matrixlist <- readmodeldata(afileinfo, yfileinfo = NULL, bmvars, bmmuttype, c("genename", bmvars , "functypecode"), qnvars, functypecodelevel = NULL, qnvarimpute=c(NA,NA), fixmusd= fixmusdfile) # note normalization of qnvar is approximately right.

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

  return(list(sigmtx = sigmtx,
              at = at,
              BMvbeta = BMvbeta,
              genesp = genesp,
              ll = ll,
              ff = ff,
              yn = yn))
}












