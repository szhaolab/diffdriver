#' Title
#'
#' @param family indicates the type of E variable
#' @param Niter is the number of replication in simulation
#' @param sgdata is the annotation data
#' @param bmrpars is the background mutation rate.
#' @param ... other parameters
#' @param Nsample, total number of samples
#' @param para para is the parameters for generating the data
#' @param bmrpars a vector, each item is bmr in log scale for one nttype
#' @param betaf0, shift of mutation rate from BMR, log scale, shared in all samples.
#' @param beta_gc, effect size for functional covariates, log scale. Right now we only have one functional covariate, that is whether the mutation is missense or loss of function mutation. beta_gc[1] indicates the shift of mutation rate for missense ( coded as 7 in functype code column in `sgdata`), beta_gc[2] indicates the shift of mutation rate for loss of function ( coded as 8 in functype code column in `sgdata`). log scale.
#' @param hotspot hotspot[1] is the probability of being hotspot for a given position. hotspot[2] is the log size
#' for hotspots.
#' @return A list composed of the p-values for 8 models
#' and the parameters used in these models.
#' @export
power_compare <- function(family, Niter, sgdata, Nsample,para,bmrpars,betaf0,beta_gc,hotseq,hmm){


  m1.pvalue <- m2.pvalue <- m3.pvalue <- m4.pvalue <-
    m5.pvalue <- m6.pvalue <- m7.pvalue <- m8.pvalue <- rep(1,Niter)
  a=c()
  for (iter in 1:Niter) {
    print(paste0("Iteration: ",  iter))
    if (family=="binary"){
      simdata <- simulate_1funcv(sgdata, bmrpars, betaf0, Nsample, beta_gc, para,hotseq,hmm)
    }else{
      if (family=="cont"){
        simdata <- simulate_2funcv(sgdata, bmrpars, betaf0, Nsample, beta_gc, para,hotseq,hmm)
      }else{
        stop("Error:only binary/cont is allowed")
      }
    }

    ssgdata=simdata$annodata
    mut <- do.call(rbind, simdata$mutlist)
    bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
    hotsize <- simdata$hotsize
    e <- simdata$pheno
    e_bisect=ifelse(e>mean(e),1,0)
    funcv <- unlist(lapply(ssgdata, "[[", "functypecode"))
    ef <- simdata$efsize
    fe1 <- c(ef$beta_gc[1], ef$beta_gc[1] + ef$beta_gc[2])[as.factor(funcv)]+hotsize
    fe2<- c(ef$avbetaf1, ef$avbetaf1 + ef$avbetaf2)[as.factor(funcv)]+hotsize
    fe3 <- rep(ef$betaf1f2, length(funcv))+hotsize
    fe4 <- rep(ef$avbetaf1f2, length(funcv))+hotsize
    mr <- bmrmtx + ef$betaf0
    if (sum(mut) ==0) {next}
    res.m1 <- mlr(mut,e)
    res.m2 <- genefisher(mut,e_bisect)
    res.m3 <- genebinom(mut,e_bisect)
    res.m4 <- genelr(mut,e_bisect)
    res.m5 <-  ddmodel(mut,e, mr, fe1)
    res.m6<- ddmodel(mut,e, mr, fe2)
    res.m7 <-  ddmodel(mut,e, mr, fe3)
    res.m8 <-  ddmodel(mut,e, mr, fe4)
    m1.pvalue[iter] <-  res.m1$pvalue
    m2.pvalue[iter] <-  res.m2$pvalue
    m3.pvalue[iter] <-  res.m3$pvalue
    m4.pvalue[iter] <-  res.m4$pvalue
    m5.pvalue[iter] <-  res.m5$pvalue
    m6.pvalue[iter] <-  res.m6$pvalue
    m7.pvalue[iter] <-  res.m7$pvalue
    m8.pvalue[iter] <-  res.m8$pvalue
    parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
    a=rbind(a,parameters)
  }
   return(list("parameters"=a, "m1.pvalue" =m1.pvalue, "m2.pvalue" =m2.pvalue,"m3.pvalue" =m3.pvalue,"m4.pvalue" =m4.pvalue,
               "m5.pvalue" =m5.pvalue,"m6.pvalue" =m6.pvalue,
               "m7.pvalue" =m7.pvalue,"m8.pvalue" =m8.pvalue,"res.m1"=res.m1, "res.m2"=res.m2, "res.m3"=res.m3, "res.m4"=res.m4,
               "res.m5"=res.m5, "res.m6"=res.m6, "res.m7"=res.m7, "res.m8"=res.m8))

  #return(list("parameters"=a, "m1.pvalue" =m1.pvalue, "m2.pvalue" =m2.pvalue,"m3.pvalue" =m3.pvalue,"m4.pvalue" =m4.pvalue))
  }
