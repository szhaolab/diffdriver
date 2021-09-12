library(data.table)

# BMR (used in simulation demo)
BMparsfile <- paste0("~/cancer_somatic/data_run/combined_20170526_5/UCS","/","UCS","_parameters_BMvar.Rdata")
load(BMparsfile)
Totalnttype <- 9
BMR <- exp(BMpars$fullpars[1:Totalnttype])/50

# BMRlist (used in diffdriver demo)
bmrdir <- "/home/simingz/cancer_somatic/data_run/combined_20170526_5/"
BMRlist <- list()
for (label in c("UCS")) {
  bmr1f <- paste0(bmrdir, label, "/",label)
  load(paste0(bmr1f,"_parameters_BMvar.Rdata"))
  load(paste0(bmr1f,"_BM_y_mu_g.Rdata"))
  BMRlist[[label]] <- list("BMpars" = BMpars, "Y_g_s_all" = Y_g_s_all, "Mu_g_s_all" = Mu_g_s_all, "nsyn" =  sum(Y_g_s_all$y)) # nsyn: total number of syn mut corresponding to this estimate
}

# Fe (used in simulation demo)
drivermapsdir <- "~/cancer_somatic/maps/"
load(paste0(drivermapsdir, "param/parmASHmean.Rdata"))
load(paste0(drivermapsdir, "param/colmu_sd_funct78.Rdata"))

Fe <- parmASHmean$TSG["functypecode8"]/allsd["functypecode8"] # effect sizes
Fe <- log(exp(Fe)*2)


# Fpars (used in diffdriver demo)
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

usethis::use_data(BMR, BMRlist, Fe, Fpars, overwrite = T)
