BMparsfile <- paste0("~/cancer_somatic/data_run/combined_20170526_5/UCS","/","UCS","_parameters_BMvar.Rdata")
load(BMparsfile)
Totalnttype <- 9
BMR <- exp(BMpars$fullpars[1:Totalnttype])/50

drivermapsdir <- "~/cancer_somatic/maps/"
load(paste0(drivermapsdir, "param/parmASHmean.Rdata"))
load(paste0(drivermapsdir, "param/colmu_sd_funct78.Rdata"))

Fe <- parmASHmean$TSG["functypecode8"]/allsd["functypecode8"] # effect sizes
Fe <- log(exp(Fe)*2)


TSGpars <- parmASHmean[[1]]
OGpars <- parmASHmean[[2]]

Fpars <- list()
OGs <- read.table(paste0(drivermapsdir,"param/OG.txt"))
TSGs <- read.table(paste0(drivermapsdir, "param/TSG.txt"))
for (g in OGs[,1]){
  Fpars[[g]] <- parmASHmean[[2]]
}
for (g in TSGs[,1]){
  Fpars[[g]] <- parmASHmean[[1]]
}

usethis::use_data(BMR, Fe, Fpars)
