<<<<<<< HEAD
# Global variable Totalnttype, Outputdir
=======
#' Global variable Totalnttype, Outputdir
>>>>>>> 94530f490d469c3e87bbf5d6ceb53c54b88b60f5
#' @import data.table
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
    ddm_j <- fread(readfile, header=TRUE, na.strings=".", colClasses= colc1, select=varlist)
  } else {
    inheader <- colnames(fread(paste("head -n 1",readfile), header=T))
    ddm_j <- fread(paste("LC_ALL=C grep -wf", genesubset, readfile), header=F, na.strings=".", colClasses= colc2,select= which(inheader %in% varlist), sep="\t")
    if (is.null(varlist)){
      colnames(ddm_j) <- inheader
    } else {colnames(ddm_j) <- inheader[inheader %in% varlist]}
  }
  return(ddm_j)
}

<<<<<<< HEAD
=======




>>>>>>> 94530f490d469c3e87bbf5d6ceb53c54b88b60f5
#' @import data.table
#' @export
ddmread <- function(afileinfo, yfileinfo, selectvars, selectmuttype, readinvars = NULL, genesubset=NULL){
  matrixlist <- list()
  selectvars <- selectvars[selectvars !="nttypecode"] # nttypecode will always be included by reading data by type
<<<<<<< HEAD
  for (j in 1:Totalnttype){
=======
  for (j in 1:totalnttype){
>>>>>>> 94530f490d469c3e87bbf5d6ceb53c54b88b60f5
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
<<<<<<< HEAD
    matrixlist[[j]] <- list(anno, y, genename)
=======
    nttypecode <- dataselect[,"nttypecode",with=F]
    matrixlist[[j]] <- list(anno, y, genename,nttypecode)
>>>>>>> 94530f490d469c3e87bbf5d6ceb53c54b88b60f5
  }
  return(matrixlist)
}

<<<<<<< HEAD
=======

>>>>>>> 94530f490d469c3e87bbf5d6ceb53c54b88b60f5
#' @import data.table
mlcoluniq <- function(matrixlist, checkcols = c("functypecode")){
  # check if all column in one dt has at least two values, if not, stop script
  # check if all column in different lists have the same set of unique values
  uniqvlist <- list()
<<<<<<< HEAD
  for (j in 1:Totalnttype){
=======
  for (j in 1:totalnttype){
>>>>>>> 94530f490d469c3e87bbf5d6ceb53c54b88b60f5
    anno <- matrixlist[[j]][[1]]
    uniq <- lapply(anno, unique)
    uniqvlist[[j]] <- lapply(uniq, sort, na.last = NA)
    uniqc <- lapply(uniq, function (x) length(x[!is.na(x)]) )
    if (length(uniq[uniqc < 2]) > 0){
      badcol <- names(uniq[uniqc < 2])
      warning(paste(c("the annotation data has only one value:", badcol, "nttype", j), collapse =" "))
    }
  }
  for (c in checkcols){
    if (length(unique(lapply(uniqvlist, `[[`, c))) > 1) stop(paste( "col", c, "have different values in different nttype"))
    print(paste(c("col",c, unlist(unique(lapply(uniqvlist, `[[`, c)))), collapse = " "))
  }
}

#' @import data.table
dtcoldup <- function(dt){
  # check if any column is duplicate, if yes, stop script
  # Note this is a rough way (check sum)
  colsum <- lapply(dt, sum , na.rm = TRUE)
  colsum <- unlist(colsum)
  if (length(unique(colsum)) < length(colsum)){
    dname <- names(colsum[colsum %in% unique(colsum[duplicated(colsum)])])
    warning(paste(c("duplicated columns found", dname), collapse = " "))
  }
}

#' @import data.table
ddmcode <- function(matrixlist, selectvars, functypecodelevel = NULL){
  # add intercept, if have functtypecode, then code and move to the front.
  selectvars <- selectvars[selectvars !="nttypecode"]
  print("coding...")
  mlcoluniq(matrixlist)
  for (j in 1:Totalnttype){
    anno <- matrixlist[[j]][[1]]
    y <- matrixlist[[j]][[2]]
    if ("functypecode" %in% selectvars){
      if (!is.null(functypecodelevel)) anno[,functypecode:= relevel(as.factor(anno[["functypecode"]]), ref = functypecodelevel)]
      newfunctype <- model.matrix( ~ functypecode, model.frame(~ functypecode, anno, na.action=na.pass ))
      anno[,functypecode:=NULL]
      newanno <- data.table(newfunctype)
    } else {
      newanno <- data.table("(Intercept)"=rep(1,dim(y)[1]))
    }
    newanno[, colnames(anno):= anno]
    dtcoldup(newanno)
    matrixlist[[j]][[1]] <- newanno
  }
  return(matrixlist)
}

<<<<<<< HEAD
=======



>>>>>>> 94530f490d469c3e87bbf5d6ceb53c54b88b60f5
#' @import data.table matrixStats
ddmprocess <- function(matrixlist, qnvars = c("expr","repl","hic"), qnvarimpute=c(-1.8,0.3), cvarimpute = 0, fixmusd=NULL){
  # for vars in qnvars, will plug in normal distributed values for missing values using a normal distribution with width of qnvarimpute[2]  and a downshift of the mean by qnvarimpute[1] compared to distribution of all.
  # for others, should all be categorical and will assign the lowest level (which is 0 under default) for missing values and then normalize.
  print("processing ...")
  print("for qnvars, filling in missing values ...")
  print("for cvars (0/1 categories), filling in missing values ...")
<<<<<<< HEAD
  for (j in 1:Totalnttype){
=======
  for (j in 1:totalnttype){
>>>>>>> 94530f490d469c3e87bbf5d6ceb53c54b88b60f5
    anno <- matrixlist[[j]][[1]]
    innames <- colnames(anno)
    inqnvars <- innames[innames %in% qnvars]
    cvars <- innames[!(innames %in% qnvars)]
    cvars <- cvars[!(cvars %in% c("(Intercept)", "chrom","start","end","ref","alt"))]
    if (length(inqnvars) > 0){
      genename <- matrixlist[[j]][[3]]
      for (qnvar in inqnvars){
        NAgns <- unique(genename[which(is.na(anno[[qnvar]]))])[["genename"]]
        vmean <- if (qnvar != "repl") qnvarimpute[1] else -qnvarimpute[1]
        v <- rnorm(length(NAgns), vmean, qnvarimpute[2]) # if positions share same genename, use same inputed value
        NAgnc <- table(genename[which(is.na(anno[[qnvar]]))])
        vall <- rep(v, NAgnc[NAgns])
        set(anno, which(is.na(anno[[qnvar]])), qnvar, vall)
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
    rm(anno, genename);gc()
  }
  if (is.null(fixmusd)){
    print("no mean or standard deviation file provided, will calculate from data ...")
    allsum <- colSums(do.call(rbind,lapply(matrixlist, function(x) x[[1]][,lapply(.SD,sum)])))
    alllen <-  colSums(do.call(rbind,lapply(matrixlist, function(x) dim(x[[1]]))))[1]
    allmu <- allsum/alllen
    allvarlist <- list()
<<<<<<< HEAD
    for (j in 1:Totalnttype){
=======
    for (j in 1:totalnttype){
>>>>>>> 94530f490d469c3e87bbf5d6ceb53c54b88b60f5
      allvarlist[[j]] <- colVars(matrixlist[[j]][[1]],center=allmu) * dim(matrixlist[[j]][[1]])[1]
      gc()
    }
    allsd <- sqrt(colSums(do.call(rbind, allvarlist))/(alllen-1))
    save(allmu, allsd, file=paste0(Outputbase,"colmu_sd_funcv.Rdata"))
  } else {
    load(fixmusd)
  }
  print("normalizing categorical variables in annotation matrix ...")
<<<<<<< HEAD
  for (j in 1:Totalnttype){
=======
  for (j in 1:totalnttype){
>>>>>>> 94530f490d469c3e87bbf5d6ceb53c54b88b60f5
    anno <- matrixlist[[j]][[1]]
    for (cvar in cvars) set(anno, j = cvar, value = (anno[[cvar]] - allmu[cvar])/allsd[cvar])
    gc()
  }
  return(matrixlist)
}

#' @import data.table
splitddm <- function(matrixlist, cpgenelist){
  # divide matrixlist into groups defined in cpgenelist, if not in any group in cpggenelist, it is put in the last group. need to optimize almost hit 20G
  outmatrixlist <- vector("list", length(cpgenelist)+1)
<<<<<<< HEAD
  for (j in 1:Totalnttype){
=======
  for (j in 1:totalnttype){
>>>>>>> 94530f490d469c3e87bbf5d6ceb53c54b88b60f5
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

#' @import data.table
#' @export
<<<<<<< HEAD
=======


>>>>>>> 94530f490d469c3e87bbf5d6ceb53c54b88b60f5
readmodeldata <- function(afileinfo, yfileinfo, selectvars, selectmuttype, readinvars = NULL,qnvars = c("expr","repl","hic"),functypecodelevel = NULL,qnvarimpute=c(-1.8,0.3), cvarimpute = 0, genesubset=NULL, fixmusd=NULL){
  # read data for a subset of genes. genesubset is a file containing gene names, one gene name per line.
  # fixmusd is a .Rd file when loaded contain allmu and allsd variables to be used in ddmprocess
  rawmlist <- ddmread(afileinfo, yfileinfo, selectvars, selectmuttype, readinvars, genesubset)
  rawmlist_code <- ddmcode(rawmlist, selectvars, functypecodelevel)
  rm(rawmlist); gc()
  matrixlist <- ddmprocess(rawmlist_code, qnvars, qnvarimpute, cvarimpute, fixmusd)
  return(matrixlist)
}

<<<<<<< HEAD
#' @import data.table
#' @export
matrixlistToGLM <- function(matrixlist, chrposmatrixlist, BMpars, mu_g_s, y_g_s, fixpars = NULL){
  GLMlist <- list()
  chrposlist <- list()
  BMcol <- names(BMpars$fullpars)[-c(1:Totalnttype,length(BMpars$fullpars)-1,length(BMpars$fullpars))]
=======





#' Compute the bmr
#'
#' @param matrixlist
#' @param chrposmatrixlist
#' @param BMpars
#' @param mu_g_s
#' @param y_g_s
#' @param fixpars
#'
#' @return
#' @export
#'
#' @examples
matrixlistToGLM_sig <- function(matrixlist, chrposmatrixlist, BMpars, mu_g_s, y_g_s, fixpars = NULL){
  GLMlist <- list()
  chrposlist <- list()
  BMcol <- names(BMpars$fullpars)[-c(1:totalnttype,length(BMpars$fullpars)-1,length(BMpars$fullpars))]
>>>>>>> 94530f490d469c3e87bbf5d6ceb53c54b88b60f5
  fixcol <- names(fixpars)
  alpha <- BMpars$fullpars["alpha"]
  lambdaPM_g_s <- data.table(agg_var = mu_g_s$agg_var, lambdaPE = log((y_g_s$y + alpha)/(mu_g_s$V1 + alpha)), key = "agg_var")
  BMvbeta <- BMpars$fullpars[1:(length(BMpars$fullpars)-2)]
<<<<<<< HEAD
  for (j in seq(1:Totalnttype)){
    BMvbetasub <- convertbeta(j, BMvbeta)
    BManno  <- matrixlist[[j]][[1]][,c("(Intercept)",BMcol), with=F]
    fixanno <- matrixlist[[j]][[1]][,fixcol, with =F]
    chrposanno <- chrposmatrixlist[[j]][[1]]
    genename <- matrixlist[[j]][[3]]
=======
  for (j in seq(1:totalnttype)){
    BMvbetasub <- convertbeta_sig(j, BMvbeta)
    #BManno  <- matrixlist[[j]][[1]][,c("(Intercept)",BMcol), with=F]
    BManno  <- matrixlist[[j]][[1]][,BMcol, with=F]
    fixanno <- matrixlist[[j]][[1]][,fixcol, with =F]
   # chrposanno <- cbind(chrposmatrixlist[[j]][[1]],chrposmatrixlist[[j]][[4]])
 chrposanno <- chrposmatrixlist[[j]][[1]]
genename <- matrixlist[[j]][[3]]
>>>>>>> 94530f490d469c3e87bbf5d6ceb53c54b88b60f5
    genename2 <- chrposmatrixlist[[j]][[3]]
    if (!identical(genename,genename2)) stop("matrixlist and chrpos not matching!")
    chrposanno[,"genename" := genename]
    lambdaPMit <- lambdaPM_g_s[genename]
    muit <- as.matrix(BManno) %*% BMvbetasub
    fixit <- as.numeric(as.matrix(fixanno) %*% as.numeric(fixpars))
    out <- matrixlist[[j]][[1]]
    out[, baseline:= rowSums(cbind(muit, lambdaPMit$lambdaPE, fixit))]
    out[, y:= matrixlist[[j]][[2]]]
<<<<<<< HEAD
    out[, c("(Intercept)",BMcol,fixcol):= NULL]
=======
    #out[, c("(Intercept)",BMcol,fixcol):= NULL]
>>>>>>> 94530f490d469c3e87bbf5d6ceb53c54b88b60f5
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

<<<<<<< HEAD
#'
convertbeta <- function(j, vbeta){
  tbeta <- vbeta[j]
  ibeta <- vbeta[-c(1:Totalnttype)]
  vbetasub <- c(tbeta,ibeta)
  return(vbetasub)
}
=======
#' @param j
#' @param vbeta
convertbeta_sig <- function(j, vbeta){
  tbeta <- vbeta[j]
  ibeta <- vbeta[-c(1:totalnttype)]
  vbetasub <- c(tbeta,ibeta)
  return(ibeta)
}




#' Functional mutation
#'
#' @param matrixlist
#' @param chrposmatrixlist
#' @param BMpars
#' @param mu_g_s
#' @param y_g_s
#' @param fixpars
#'
#' @return A list
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


#' @param j
#' @param vbeta
convertbeta <- function(j, vbeta){
  tbeta <- vbeta[j]
  ibeta <- vbeta[-c(1:totalnttype)]
  vbetasub <- c(tbeta,ibeta)
  return(vbetasub)
}





>>>>>>> 94530f490d469c3e87bbf5d6ceb53c54b88b60f5
