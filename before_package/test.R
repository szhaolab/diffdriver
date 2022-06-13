## cmodel.frac and ddmodel should give same estimation results
source('~/cancer_pheno/cancer_pheno/ddmodel.R')
source('~/cancer_pheno/cancer_pheno/countLRT.R')
source('~/cancer_pheno/cancer_pheno/simulate.R')

load(paste0("/home/simingz/cancer_pheno/data_run/simulation_2019-04-29/sgdata.Rd"))
load("~/cancer_somatic/data_run/combined_20170526_5/UCS/UCS_parameters_BMvar.Rdata")
bmrpars <- BMpars$fullpars[1:9] - log(50)
simdata <- simulate_1funcv(sgdata, bmrpars, betaf0=0.5, Nsample=1000, Nc=200, beta_gc=c(0,1), fracc=0.8, fracn=0.2)

BMR <- exp(simdata$bmrpars)
cmodel.frac(simdata$mutbyt, simdata$pheno, simdata$annodata, c(simdata$efsize$avbetaf0, simdata$efsize$avbetaf1)) # diffDriver
cmodel.frac(simdata$mutbyt, simdata$pheno, simdata$annodata, c(simdata$efsize$avbetaf0f1, 0)) # diffDriver-baseline

mutmtx <- do.call(rbind, simdata$mutbyt)
e <- simdata$pheno
bmrmtx <- do.call(rbind, lapply(1:length(simdata$mutbyt), function(t) matrix(simdata$bmrpars[t], nrow = dim(simdata$mutbyt[[t]])[1], ncol = dim(simdata$mutbyt[[t]])[2])))
ft <- do.call(rbind, simdata$annodata)$functypecode
fe1 <- c("7" = simdata$efsize$avbetaf0, "8" = simdata$efsize$avbetaf0 + simdata$efsize$avbetaf1)[ft]
fe2 <- rep(simdata$efsize$avbetaf0f1, dim(mutmtx)[1])

ddmodel(mutmtx, e, bmrmtx, fe1)  # diffDriver
ddmodel(mutmtx, e, bmrmtx, fe2) # diffDriver-baseline