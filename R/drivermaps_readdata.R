#' Global variable Totalnttype, Outputdir
#'
#' @param fileinfo Anno file information
#' @param j the number of nytype
#' @param varlist Variables of interest
#' @param genesubset Genes of interest
#'
#' @import data.table
#' @noRd
ddmread_j <-function(fileinfo, j, varlist = NULL, genesubset = NULL){
  # genesubset is a file containing genenames, each line has one gene name
  # working in data.table > 1.10
  readfile <- gsub("XXX", j, fileinfo[["file"]])
  if (length(unique(varlist)) == length(fileinfo[["header"]])){varlist = NULL}
  headertemp <- fileinfo[["header"]][fileinfo[["header"]] %in% varlist]
  coltypetemp <- fileinfo[["coltype"]][fileinfo[["header"]] %in% varlist]
  colc1 <- lapply(unique(coltypetemp),function(x) headertemp[which(coltypetemp==x)])
  names(colc1) <- unique(coltypetemp)

  headertemp2 <- paste0("V", 1:length(fileinfo[["header"]]))[fileinfo[["header"]] %in% varlist]
  colc2 <- lapply(unique(coltypetemp),function(x) headertemp2[which(coltypetemp==x)])
  names(colc2) <- unique(coltypetemp)

  if (is.null(genesubset)) {
    ddm_j <- suppressMessages(data.table::fread(readfile, header=TRUE, na.strings=".", colClasses= colc1, select=varlist))
  } else {
    inheader <- suppressMessages(colnames(data.table::fread(paste("head -n 1",readfile), header=T)))
    ddm_j <- suppressMessages(data.table::fread(paste("LC_ALL=C grep -wf", genesubset, readfile), header=F, na.strings=".", colClasses= colc2,select= which(inheader %in% varlist), sep="\t"))
    if (is.null(varlist)){
      colnames(ddm_j) <- inheader
    } else {colnames(ddm_j) <- inheader[inheader %in% varlist]}
  }
  return(ddm_j)
}



ddmread <- function(afileinfo, yfileinfo, selectvars, selectmuttype, readinvars = NULL, genesubset=NULL){
  matrixlist <- list()
  selectvars <- selectvars[selectvars !="nttypecode"] # nttypecode will always be included by reading data by type
  for (j in 1:totalnttype){
    dataall <- ddmread_j(afileinfo, j, varlist = readinvars, genesubset=genesubset)
    if (is.null(yfileinfo)){
      # if no Yfileinfo use 0 for all positions
      dataall[,y:= 0]
    } else {
      yall <- ddmread_j(yfileinfo, j, varlist = "y", genesubset=genesubset)
      suppressWarnings(dataall[,y:=yall]) #  Invalid .internal.selfref warning
    }
    dataselect <- dataall[eval(parse(text =selectmuttype))]
    anno <- dataselect[ , selectvars, with = F]
    y <- dataselect[ , "y", with = F]
    genename <- dataselect[ , "genename", with = F]
    nttypecode <- dataselect[,"nttypecode",with=F]
    matrixlist[[j]] <- list(anno, y, genename,nttypecode)
  }
  return(matrixlist)
}


#'  Add intercept, if have functtypecode, then code and move to the front.
#'  different from driverMAPS! only allows
#'  functtypecode =7 ||8 when functypecode is included in selectvars.
ddmcode <- function(matrixlist, selectvars, functypecodelevel = NULL){
  selectvars <- selectvars[selectvars !="nttypecode"]
  print("coding...")
  for (j in 1:totalnttype){
    anno <- matrixlist[[j]][[1]]
    y <- matrixlist[[j]][[2]]
    newanno <- data.table("(Intercept)"=rep(1,dim(y)[1]))
    if ("functypecode" %in% selectvars){
      newanno[anno$functypecode ==7, functypecode8:= 0]
      newanno[anno$functypecode ==8, functypecode8:= 1]
      anno[,functypecode:=NULL]
    }
    newanno[, colnames(anno):= anno]
    matrixlist[[j]][[1]] <- newanno
  }
  return(matrixlist)
}


ddmprocess <- function(matrixlist, qnvars = c("expr","repl","hic"), qnvarimpute=c(0,0), cvarimpute = 0, fixmusd=NULL){
print("processing ...")
  print("for qnvars, filling in missing values ...")
  print("for cvars (0/1 categories), filling in missing values ...")
  for (j in 1:totalnttype){
    anno <- matrixlist[[j]][[1]]
    innames <- colnames(anno)
    inqnvars <- innames[innames %in% qnvars]
    cvars <- innames[!(innames %in% qnvars)]
    cvars <- cvars[!(cvars %in% c("(Intercept)", "chrom","start","end","ref","alt"))]
    if (length(inqnvars) > 0){
      # a little different from driverMAPS
      for (qnvar in inqnvars){
        set(anno, which(is.na(anno[[qnvar]])), qnvar, qnvarimpute[1])
      }
    }
    # for cvars (0/1 categories), filling in missing values ...
    if (length(cvars) > 0) {
      for (cvar in cvars){
        set(anno, which(is.na(anno[[cvar]])), j = cvar, value = cvarimpute)
      }
    }
    matrixlist[[j]][[1]] <- anno[complete.cases(anno),]
    matrixlist[[j]][[2]] <- matrixlist[[j]][[2]][complete.cases(anno),]
    matrixlist[[j]][[3]] <- matrixlist[[j]][[3]][complete.cases(anno),]
    rm(anno);gc()
  }
  if (is.null(fixmusd)){
    print("no mean or standard deviation file provided, will calculate from data ...")
    allsum <- colSums(do.call(rbind,lapply(matrixlist, function(x) x[[1]][,lapply(.SD,sum)])))
    alllen <-  colSums(do.call(rbind,lapply(matrixlist, function(x) dim(x[[1]]))))[1]
    allmu <- allsum/alllen
    allvarlist <- list()
    for (j in 1:totalnttype){
      # Note this does't work with matrixStats anymore
      allvarlist[[j]] <- matrixStats::colVars(matrixlist[[j]][[1]],center=allmu) * dim(matrixlist[[j]][[1]])[1]
      gc()
    }
    allsd <- sqrt(colSums(do.call(rbind, allvarlist))/(alllen-1))
    save(allmu, allsd, file=paste0(outputbase,"colmu_sd_funcv.Rdata"))
  } else {
    load(fixmusd)
  }
  print("normalizing categorical variables in annotation matrix ...")
  for (j in 1:totalnttype){
    anno <- matrixlist[[j]][[1]]
    for (cvar in cvars) set(anno, j = cvar, value = (anno[[cvar]] - allmu[cvar])/allsd[cvar])
    gc()
  }
  return(matrixlist)
}


#' divide matrixlist into groups defined in cpgenelist,
#' if not in any group in cpggenelist,
#' it is put in the last group. need to optimize almost hit 20G
#'
#' @param matrixlist
#' @param cpgenelist
#' @return A list
#' @noRd
splitddm <- function(matrixlist, cpgenelist){
  outmatrixlist <- vector("list", length(cpgenelist)+1)
  for (j in 1:totalnttype){
    annoall <- matrixlist[[j]][[1]]
    yall <- matrixlist[[j]][[2]]
    geneall <- matrixlist[[j]][[3]]
    for (m in 1:length(outmatrixlist)) {
      if (m != length(outmatrixlist)) {
        anno <- annoall[geneall[["genename"]] %in% cpgenelist[[m]]]
        y <- yall[geneall[["genename"]] %in% cpgenelist[[m]]]
        gene <- geneall[genename %in% cpgenelist[[m]]]
      } else {
        anno <- annoall[!(geneall[["genename"]] %in% as.matrix(unlist(cpgenelist)))]
        y <- yall[!(geneall[["genename"]] %in% as.matrix(unlist(cpgenelist)))]
        gene <- geneall[!(genename %in% as.matrix(unlist(cpgenelist)))]
      }
      outmatrixlist[[m]][[j]] <- list(anno, y,gene)
      gc()
    }
  }
  return(outmatrixlist)
}


readmodeldata <- function(afileinfo, yfileinfo, selectvars, selectmuttype, readinvars = NULL,qnvars = c("expr","repl","hic"),functypecodelevel = NULL,qnvarimpute=c(-1.8,0.3), cvarimpute = 0, genesubset=NULL, fixmusd=NULL){
  # read data for a subset of genes. genesubset is a file containing gene names, one gene name per line.
  # fixmusd is a .Rd file when loaded contain allmu and allsd variables to be used in ddmprocess
  rawmlist <- ddmread(afileinfo, yfileinfo, selectvars, selectmuttype, readinvars, genesubset)
  rawmlist_code <- ddmcode(rawmlist, selectvars, functypecodelevel)
  rm(rawmlist); gc()
  matrixlist <- ddmprocess(rawmlist_code, qnvars, qnvarimpute, cvarimpute, fixmusd)
  return(matrixlist)
}



matrixlistToGLM <- function(matrixlist, chrposmatrixlist, BMpars, mu_g_s, y_g_s, fixpars = NULL){
  GLMlist <- list()
  chrposlist <- list()
  BMcol <- names(BMpars$fullpars)[-c(1:totalnttype,length(BMpars$fullpars)-1,length(BMpars$fullpars))]
  fixcol <- names(fixpars)
  alpha <- BMpars$fullpars["alpha"]
  lambdaPM_g_s <- data.table(agg_var = mu_g_s$agg_var, lambdaPE = log((y_g_s$y + alpha)/(mu_g_s$V1 + alpha)), key = "agg_var")
  BMvbeta <- BMpars$fullpars[1:(length(BMpars$fullpars)-2)]
  for (j in seq(1:totalnttype)){
    BMvbetasub <- convertbeta(j, BMvbeta)
    BManno  <- matrixlist[[j]][[1]][,c("(Intercept)",BMcol), with=F]
    fixanno <- matrixlist[[j]][[1]][,fixcol, with =F]
    chrposanno <- cbind(data.table(chrposmatrixlist[[j]][[1]]),data.table(chrposmatrixlist[[j]][[4]]))
# chrposanno <- chrposmatrixlist[[j]][[1]]
    genename <- matrixlist[[j]][[3]]
    genename2 <- chrposmatrixlist[[j]][[3]]
    if (!identical(genename,genename2)) stop("matrixlist and chrpos not matching!")
    chrposanno[,"genename" := genename]
    lambdaPMit <- lambdaPM_g_s[genename]
    muit <- as.matrix(BManno) %*% BMvbetasub
    fixit <- as.numeric(as.matrix(fixanno) %*% as.numeric(fixpars))
    out <- matrixlist[[j]][[1]]
    out[, baseline:= rowSums(cbind(muit, lambdaPMit$lambdaPE, fixit))]
    out[, y:= matrixlist[[j]][[2]]]
    #out[, c("(Intercept)",BMcol,fixcol):= NULL]
    GLMlist[[j]] <- out
    chrposlist[[j]] <- chrposanno
  }
  glmdt <- do.call(rbind,GLMlist)
  chrposdt <- do.call(rbind,chrposlist)
  selectcol <- copy(colnames(glmdt))
  glmdt[,colnames(chrposdt):= chrposdt]
  glmdtsorted <- glmdt[order(chrom,genename,start)]
  glmdt <- glmdtsorted[,selectcol,with = F]
  glmdtchrpos <- glmdtsorted[,c("chrom","genename","start", colnames(chrposdt)), with =F]
  return(list(glmdt,glmdtchrpos))
}

convertbeta <- function(j, vbeta){
  tbeta <- vbeta[j]
  ibeta <- vbeta[-c(1:totalnttype)]
  vbetasub <- c(tbeta,ibeta)
  return(vbetasub)
}


### ----- clean up---------
matrixlistToGLM_sig <- function(matrixlist, chrposmatrixlist, BMpars, mu_g_s, y_g_s, fixpars = NULL){
  GLMlist <- list()
  chrposlist <- list()
  BMcol <- names(BMpars$fullpars)[-c(1:totalnttype,length(BMpars$fullpars)-1,length(BMpars$fullpars))]
  fixcol <- names(fixpars)
  alpha <- BMpars$fullpars["alpha"]
  lambdaPM_g_s <- data.table(agg_var = mu_g_s$agg_var, lambdaPE = log((y_g_s$y + alpha)/(mu_g_s$V1 + alpha)), key = "agg_var")
  BMvbeta <- BMpars$fullpars[1:(length(BMpars$fullpars)-2)]
  for (j in seq(1:totalnttype)){
    BMvbetasub <- convertbeta_sig(j, BMvbeta)
    #BManno  <- matrixlist[[j]][[1]][,c("(Intercept)",BMcol), with=F]
    BManno  <- matrixlist[[j]][[1]][,BMcol, with=F]
    fixanno <- matrixlist[[j]][[1]][,fixcol, with =F]
    # chrposanno <- cbind(chrposmatrixlist[[j]][[1]],chrposmatrixlist[[j]][[4]])
    chrposanno <- chrposmatrixlist[[j]][[1]]
    genename <- matrixlist[[j]][[3]]
    genename2 <- chrposmatrixlist[[j]][[3]]
    if (!identical(genename,genename2)) stop("matrixlist and chrpos not matching!")
    chrposanno[,"genename" := genename]
    lambdaPMit <- lambdaPM_g_s[genename]
    muit <- as.matrix(BManno) %*% BMvbetasub
    fixit <- as.numeric(as.matrix(fixanno) %*% as.numeric(fixpars))
    out <- matrixlist[[j]][[1]]
    out[, baseline:= rowSums(cbind(muit, lambdaPMit$lambdaPE, fixit))]
    out[, y:= matrixlist[[j]][[2]]]
    #out[, c("(Intercept)",BMcol,fixcol):= NULL]
    GLMlist[[j]] <- out
    chrposlist[[j]] <- chrposanno
  }
  glmdt <- do.call(rbind,GLMlist)
  chrposdt <- do.call(rbind,chrposlist)
  selectcol <- copy(colnames(glmdt))
  glmdt[,colnames(chrposdt):= chrposdt]
  glmdtsorted <- glmdt[order(chrom,genename,start)]
  glmdt <- glmdtsorted[,selectcol,with = F]
  glmdtchrpos <- glmdtsorted[,c("chrom","genename","start", colnames(chrposdt)), with =F]
  return(list(glmdt,glmdtchrpos))
}

convertbeta_sig <- function(j, vbeta){
  tbeta <- vbeta[j]
  ibeta <- vbeta[-c(1:totalnttype)]
  vbetasub <- c(tbeta,ibeta)
  return(ibeta)
}





