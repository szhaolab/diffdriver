#' @title Simulation function
#'
#' @param family indicates the type of E variable
#' @param Niter is the number of replication in simulation
#' @param sgdata is the annotation data
#' @param bmrpars is the background mutation rate.
#' @param betaf0 gene effect mutation
#' @param Nsample sample size of simulation
#' @param beta_gc the effects of missense and nonsense
#' @param par parameters for the simuation
#'
#' @return A list composed of the p-values for 8 models
#'        and the parameters used in these models.
#' @export
power_compare <- function(family="binary", Niter=200, sgdata, bmrpars, betaf0=0,Nsample=1000,beta_gc=c(0,1.2),par){
  m1.pvalue <- m2.pvalue <- m3.pvalue <- m4.pvalue <-
    m5.pvalue <- m6.pvalue <- m7.pvalue <- m8.pvalue <-rep(1,Niter)
  a=c()
  for (iter in 1:Niter) {
    print(paste0("Iteration: ",  iter))
    if (family=="binary"){
      simdata <- simulate_1funcv(sgdata, bmrpars, betaf0, Nsample, beta_gc, par)
    }else{
    if (family=="cont"){
      simdata <- simulate_2funcv(sgdata, bmrpars, betaf0, Nsample, beta_gc, par)
    }else{
      stop("Error:only binary/cont is allowed")
    }
    }
    ssgdata=simdata$annodata
    mut <- do.call(rbind, simdata$mutlist)
    bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
    e <- simdata$pheno
    e_bisect=ifelse(e>mean(e),1,0)
    if (sum(mut) ==0) {next}
    res.m1 <- mlr(mut,e)
    res.m2 <- genefisher(mut,e_bisect)
    res.m3 <- genebinom(mut,e_bisect)
    res.m4 <- genelr(mut,e_bisect)
    funcv <- unlist(lapply(ssgdata, "[[", "functypecode"))
    ef <- simdata$efsize
    fe1 <- c(ef$beta_gc[1], ef$beta_gc[1] + ef$beta_gc[2])[as.factor(funcv)]
    fe2<- c(ef$avbetaf1, ef$avbetaf1 + ef$avbetaf2)[as.factor(funcv)]
    fe3 <- rep(ef$betaf1f2, length(funcv))
    fe4 <- rep(ef$avbetaf1f2, length(funcv))
    mr <- bmrmtx + ef$betaf0
    res.m5 <-  ddmodel(mut,e, mr, fe1)
    res.m6<- ddmodel(mut,e, mr, fe2)
    parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
    a=rbind(a,parameters)
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
  }
  return(list("parameters"=a, "m1.pvalue" =m1.pvalue, "m2.pvalue" =m2.pvalue,"m3.pvalue" =m3.pvalue,"m4.pvalue" =m4.pvalue,
              "m5.pvalue" =m5.pvalue,"m6.pvalue" =m6.pvalue,
              "m7.pvalue" =m7.pvalue,"m8.pvalue" =m8.pvalue))
}
