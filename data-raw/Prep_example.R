library(data.table)

# BMR (used in simulation demo)
BMparsfile <- paste0("~/cancer_somatic/data_run/combined_20170526_5/UCS","/","UCS","_parameters_BMvar.Rdata")
load(BMparsfile)
Totalnttype <- 9
BMR <- exp(BMpars$fullpars[1:Totalnttype])/50


# Fe (used in simulation demo)
drivermapsdir <- "~/cancer_somatic/maps/"
load(paste0(drivermapsdir, "param/parmASHmean.Rdata"))
load(paste0(drivermapsdir, "param/colmu_sd_funct78.Rdata"))

Fe <- parmASHmean$TSG["functypecode8"]/allsd["functypecode8"] # effect sizes
Fe <- log(exp(Fe)*2)

# Data for one single gene

library(diffdriver)
drivermapsdir <- "~/cancer_somatic/maps/"

sg <- "ERBB3"
Totalnttype <- 9
Adirbase < paste0(drivermapsdir, "quicktest_data/")
Afileinfo <- list(file = paste(Adirbase, "nttypeXXX_annodata.txt", sep=""),
                  header = c("chrom","start","end","ref","alt","genename","functypecode","nttypecode","expr","repl","hic","mycons","sift","phylop100","MA","ssp","wggerp"),
                  coltype = c("character","numeric","numeric","character","character","character","character","factor","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))

dataall <- list()
sgdata <- list()
for (j in 1:Totalnttype){
  dataall[[j]] <- ddmread_j(Afileinfo, j, varlist = c("chrom","start","genename","functypecode","nttypecode"))
  sgdata[[j]] <- dataall[[j]][(functypecode==7 | functypecode==8)& genename == sg]
}
print(sgdata[[1]])

usethis::use_data(BMR, BMRlist, Fe, overwrite = T)
#usethis::use_data(sgdata, overwrite = T)

# Files
# Mutatiion file: from TCGA-UCS
# Phenotype: smokingcessation PRS for UCS from Dr. Jian Yang'g group.
# Hotspot: any consecutive positions with mutations. TODO: to be removed when diffdriver can do this.

