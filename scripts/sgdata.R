
annoAllg <- function(genef, mutf, phenof,j,hotf, drivermapsdir, outputdir =".", outputname = "diffdriver_results"){
  # ------- read position level information (same as in drivermaps) ----------
  adirbase <-drivermapsdir
  afileinfo <- list(file = paste(adirbase, "/nttypeXXX_annodata.txt", sep=""),
                    header = c("chrom","start","end","ref","alt","genename","functypecode","nttypecode","expr","repl","hic","mycons","sift","phylop100","MA","ssp","wggerp"),
                    coltype = c("character","numeric","numeric","character","character","character","character","factor","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))
  totalnttype <<- 9
  Totalnttype <<- 9
 bmvars <- c("nttypecode", "expr", "repl", "hic")
  bmmuttype <- "functypecode == 6 & ssp == 0"
  funcvars <- c("functypecode", "mycons", "sift", "phylop100", "MA")
  functypecodelevel <- "7"
  funcvmuttype <- "functypecode == 7 | functypecode == 8"
  readinvars <- c("genename", "ssp",bmvars, funcvars) # This is optional, if not given, then will read all columns
  qnvars = c("expr","repl","hic") # all the rest will be normalized, except for nttypecode
  outputbase <<- paste0(outputdir, outputname)
  paramdir <- paste0(drivermapsdir, "param/")
  fixmusdfile <-  paste0(paramdir, "colmu_sd_funct78.Rdata")

  allg <- read.table(genef, stringsAsFactors = F)[,1]
matrixlist <- readmodeldata(afileinfo, yfileinfo = NULL, c(bmvars,funcvars), funcvmuttype, readinvars , qnvars, functypecodelevel,qnvarimpute=c(0,0), cvarimpute = 0, genesubset=genef, fixmusd=fixmusdfile)
chrposmatrixlist <- ddmread(afileinfo, yfileinfo = NULL, c("chrom", "start","ref","alt"), funcvmuttype, c("genename", "chrom", "start", "ref", "alt", "functypecode", "ssp", "nttypecode"), genesubset=genef)

#save(matrixlist,file="matrixlist9.Rd")
#save(chrposmatrixlist,file="chrposmatrixlist9.Rd")
#load("matrixlist9.Rd")
#load("chrposmatrixlist9.Rd")
BMRlist=BMRlist$UCS

for (t in 1:length(matrixlist)){
    b1 <- which(sapply(matrixlist[[t]][[3]], grepl, pattern="[;,|]"))
    matrixlist[[t]][[3]][b1,] <- unlist(lapply(sapply(matrixlist[[t]][[3]][b1,], strsplit, split="[;,|]"), function(x) intersect(x,allg)[1]))
    b2 <- which(sapply(chrposmatrixlist[[t]][[3]], grepl, pattern="[;,|]"))
    chrposmatrixlist[[t]][[3]][b2,] <- unlist(lapply(sapply(matrixlist[[t]][[3]][b2,], strsplit, split="[;,|]"), function(x) intersect(x,allg)[1]))

    }
  #------------------------------------------------------------------------------

  # background mutation rate (bmrdt): data.table, each column correponds to one BMR label
  bmrdt <- data.table()


    matrixlisttemp <- copy(matrixlist)
    chrposmatrixlisttemp <- copy(chrposmatrixlist)
    y_g_s <- BMRlist$Y_g_s_all[allg][,1:2, with=F]
    y_g_s[is.na(y_g_s)] <- 0
    mu_g_s <- BMRlist$Mu_g_s_all[allg][,1:2, with=F]
    mu_g_s[is.na(mu_g_s)] <- 0
    glmdtall <- matrixlistToGLM(matrixlisttemp, chrposmatrixlisttemp, BMRlist$BMpars, mu_g_s, y_g_s, fixpars=NULL)
    #glmdtall <- matrixlistToGLM(matrixlist, chrposmatrixlist, BMRlist[[label]]$BMpars, mu_g_s, y_g_s, fixpars=NULL)
    bmrdt[,eval(type):= glmdtall[[1]]$baseline]
    rm(matrixlisttemp,chrposmatrixlisttemp); gc()


  # functional annotation (fanno):data.table, each row has functional annotation for each possible mutation
  fanno <- glmdtall[[1]]

  # row index (ri): chr pos ref alt
  ri <- glmdtall[[2]]
return (list=c("bmr"=bmr,"fannomatrix"=fanno,"fanno"=ri))
}
