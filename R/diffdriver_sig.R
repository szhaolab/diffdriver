# TODO: mutation hotspot, quantitative phenotype
# Notes: syn data prepared from direct join of mutation list to annodata is different from used in driverMAPS BMR estimation: 1. 5% syn are ssp are included. 2. BMR estimation from old driverMAPS run has different annodata. 3. genes without expr/hic/rep were removed in BMR driverMAPS estimation.
#' @title Run diffDriver given input files
#' @description This is the function to run diffDriver. We first need to set up: run driverMAPS for groups with potential different BMR, assign BMR labels for each sample. then BMR for each sample will be scaled based on driverMAPS results.
#' @param genef file for name of genes to be included in the analysis
#' @param mutf mutation list file, use the driverMAPS mutation input format
#' @param phenof phenptype file, SampleID <tab> Phenotype <tab> Nsyn. nsyn is number of syn mutations in this sample.
#' @param j The index of phenotype
#' @import Matrix data.table
#' @export
diffdriver_sig= function(genef, mutf, phenof,j, hotf, drivermapsdir,k=k, outputdir =".", outputname = "diffdriver_results"){
  # ------- read position level information (same as in drivermaps) ----------
  adirbase <-drivermapsdir
  afileinfo <- list(file = paste(adirbase, "TCGA-UCS_nttypeXXX_annodata.txt", sep=""),
                    header = c("chrom","start","end","ref","alt","genename","functypecode","nttypecode","expr","repl","hic","mycons","sift","phylop100","MA","ssp","wggerp"),
                    coltype = c("character","numeric","numeric","character","character","character","character","factor","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))
  totalnttype <<- 96
  Totalnttype <<- 96
  bmvars <- c("nttypecode", "expr", "repl", "hic")
  bmmuttype <- "functypecode == 6"
  funcvars <- c("functypecode", "mycons", "sift", "phylop100", "MA")
  functypecodelevel <- "7"
  funcvmuttype <- "functypecode == 7 | functypecode == 8"
  readinvars <- c("genename", "ssp", bmvars, funcvars) # This is optional, if not given, then will read all columns
  bmreadinvars <- c("genename", "ssp","nttypecode",bmvars, funcvars)
  qnvars = c("expr","repl","hic") # all the rest will be normalized, except for nttypecode
  outputbase <<- paste0(outputdir, "/", outputname)
  paramdir <- paste0(drivermapsdir, "/param/")
  fixmusdfile <-  paste0(paramdir, "colmu_sd_funct78.Rdata")

  allg <- read.table(genef, stringsAsFactors = F)[,1]
 matrixlist <- readmodeldata(afileinfo, yfileinfo = NULL, c(bmvars,funcvars), funcvmuttype, readinvars , qnvars, functypecodelevel,qnvarimpute=c(0,0), cvarimpute = 0, genesubset=genef, fixmusd=fixmusdfile)
 chrposmatrixlist <- ddmread(afileinfo, yfileinfo = NULL, c("chrom", "start","ref","alt","nttypecode"), funcvmuttype, c("genename", "chrom", "start", "ref", "alt", "functypecode", "ssp", "nttypecode"), genesubset=genef)
#save(matrixlist,file="matrixlist96.Rd")
#save(chrposmatrixlist,file="chrposmatrixlist96.Rd")
#load("matrixlist96.Rd")
#load("chrposmatrixlist96.Rd")
 BMRlist=BMRlist$UCS


  for (t in 1:length(matrixlist)){
    b1 <- which(sapply(matrixlist[[t]][[3]], grepl, pattern="[;,|]"))
    matrixlist[[t]][[3]][b1,] <- unlist(lapply(sapply(matrixlist[[t]][[3]][b1,], strsplit, split="[;,|]"), function(x) intersect(x,allg)[1]))
    b2 <- which(sapply(chrposmatrixlist[[t]][[3]], grepl, pattern="[;,|]"))
    chrposmatrixlist[[t]][[3]][b2,] <- unlist(lapply(sapply(matrixlist[[t]][[3]][b2,], strsplit, split="[;,|]"), function(x) intersect(x,allg)[1]))

  }
  #------------------------------------------------------------------------------

  # background mutation rate (bmrdt): data.table, each column correponds to one BMR label
#target=vector("list",length(BMRlist))
#names(target)<- names(BMRlist)

    matrixlisttemp <- copy(matrixlist)
    chrposmatrixlisttemp <- copy(chrposmatrixlist)
    y_g_s <- BMRlist$Y_g_s_all[allg][,1:2, with=F]
    y_g_s[is.na(y_g_s)] <- 0
    mu_g_s <- BMRlist$Mu_g_s_all[allg][,1:2, with=F]
    mu_g_s[is.na(mu_g_s)] <- 0
    glmdtall <- matrixlistToGLM(matrixlisttemp, chrposmatrixlisttemp, BMRlist$BMpars, mu_g_s, y_g_s, fixpars=NULL)
    #bmrdt[,eval(label):= glmdtall[[1]]$baseline]
    #target=cbind(glmdtall[[2]][,.(chrom,genename,start,nttypecode)],glmdtall[[1]][,.(expr,repl,hic)])
    rm(matrixlisttemp,chrposmatrixlisttemp,matrixlist,chrposmatrixlist); gc()




    bmrsig <- matrixlistToBMR(adirbase, mutf,BMRlist,k=k)




  # functional annotation (fanno):data.table, each row has functional annotation for each possible mutation
  fanno <- glmdtall[[1]]

  # row index (ri): chr pos ref alt
  ri <- glmdtall[[2]][,.(chrom,genename,start,ref,alt,nttypecode)]
  lambda_pe=bmrsig$lambda
  ll=bmrsig$ll
  ff=bmrsig$ff
  save(ll,file=paste0(outputbase,"_loadings.Rd"))
  save(ff,file=paste0(outputbase,"_factors.Rd"))
ri=plyr::join(ri,lambda_pe)
index=unique(which(is.na(ri),arr.ind = T)[,1])
fanno=fanno[-index,]
ri=ri[-index,]

# mutations (muts): data.table, with columns Chromosome, Position, Ref, Alt, SampleID
muts0 <- fread(mutf, header = T)
if (!grepl('chr', muts0$Chromosome[1], fixed = T)) {muts0$Chromosome <- paste0("chr",muts0$Chromosome)}
  # sample annotation (canno):data.table, with columns BMR label, No. syn and phenotype.
  canno0 <- as.data.table(fread(phenof, header = "auto"))
shared=intersect(muts0$SampleID,canno0$SampleID)
index1=which(canno0$SampleID %in% shared)
index2=which(muts0$SampleID %in% shared)
muts<- muts0[index2,]
canno = canno0[index1,]

  # column index (ci): sampleID
  ci <- canno[,"SampleID"]
  ci[,"cidx" := 1:dim(canno)[1]]

  # for (l in names(BMRlist)) {
  #   BMRlist[[l]][["nsyn"]] <- sum(canno$Nsyn[canno$BMRlabel == l])
  # }


##estimate bmr mu_ij
mutation=data.table()
coe=BMRlist$BMpars$fullpars[c("expr","repl","hic")]
numerator=exp(fanno$expr*coe[1]+fanno$repl*coe[2]+fanno$hic*coe[3])
for (i in 1:length(shared)) {
  print(paste("bmr for sample",i,sep=" "))
  mui=(ri$lambda)*numerator*(bmrsig$sampelsig[i,ri$nttypecode])/bmrsig$normconstant
  mutation=cbind(mutation,mui)
    }
if (any(is.na(mutation))) {stop("bmr missing")}
  riallg <- split(ri,ri$genename)
  fannoallg <- split(fanno,ri$genename)
bmrmtx <- split(mutation,ri$genename)
rm(ri,mutation,fanno)
  # run diffdriver for each gene
  res <- list()
  ## hotspot
  hotspots=read.table(file = hotf)
  hmm=readRDS(paste0(dirvermapdir, "hmmOGpar_ASHmean.rds"))
  for (g in names(bmrmtx)) {
    print(paste0("Start to process gene: ", g))
    rig <- riallg[[g]]
    rig$ridx <- 1:dim(rig)[1]
    muti <- na.omit(ci[rig[muts, on = c("chrom"= "Chromosome", "start" = "Position",  "ref" = "Ref",  "alt"= "Alt")], on = "SampleID"])
    mutmtx <- sparseMatrix(i = muti$ridx, j = muti$cidx, dims = c(max(rig$ridx), max(ci$cidx)))

    hotg= na.omit(rig[hotspots,on=c("chrom"="chrom","start"="start")])
    hotmat=rep(0,nrow(rig))
    hotmat[hotg$ridx]=1
    if (sum(mutmtx) ==0) {
      next
      }


    betaf <- Fpars[[g]][names(Fpars[[g]]) != "beta_f0"]
    betaf0 <- Fpars[[g]]["beta_f0"]
    if (is.null(betaf)){
      betaf <-  Fpars[["TP53"]][names(Fpars[["TP53"]]) != "beta_f0"]
      betaf0 <-  Fpars[["TP53"]]["beta_f0"]
    } # if OG/TSG unknown, use TSG parameters.
    ganno <- fannoallg[[g]]
    fe <- as.matrix(ganno[ ,names(betaf), with =F]) %*% betaf+ hotmat*hmm[8] + betaf0

    resg <- list()
    e=canno[[j]]
    resg <- ddmodel(mutmtx, e, bmrmtx[[g]], fe[,1])
    resg[["mlr"]] <- mlr(mutmtx, e)
    resg[["mlr.v2"]] <- mlr.v2(mutmtx, e, canno$Nsyn)
    resg[["fisher"]] <- genefisher(mutmtx, e)
    #resg[["binom"]] <- genebinom(mutmtx, e)
    #resg[["lr"]] <- genelr(mutmtx, e)
    res[[g]] <- resg

    save(canno,fe, ganno, betaf, betaf0, resg, file=paste0(paste0(outputbase,".", g, ".Rd")))
    #setEPS()
    #postscript(file=paste0(outputbase,".", g, "mut_status.eps"), width=9, height=4)
    #plot_mut(mutmtx, canno, bmrmtx, ganno)
    #dev.off()
  }

  return(res)
}
