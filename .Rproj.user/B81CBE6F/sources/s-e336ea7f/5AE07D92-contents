# readin data for one single gene

library(diffdriver)
drivermapsdir <- "~/cancer_somatic/maps/"
load_drivermaps_func(drivermapsdir)

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

usethis::use_data(sgdata, overwrite = T)
