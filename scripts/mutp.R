

#' Title
#'
#' @param gene is the gene of interest
#' @param funcanno is the functional annotations that are used in the model
#' @param Afiledir is the path to the annotation files
#' @param tau is the tuning parameter used to test the mutation to the parameter setting.
#'
#' @return
#' @export
#'
#' @examples
#mutp <- function(gene,funcanno,Afiledir,betaf0,  beta_gc){
#          mutp=vector("list",9)
#          genedata=vector("list",9)
#		for (t in 1:9) {
#			genedata[[t]]=fread(paste0(Afiledir,"/nttype",t,"_annorate.txt"))[genename==gene & functypecode==7 & functypecode==8,]
#			bmr=genedata[[t]][,nttypecode]*genedata[[t]][,geffect]*betaf0
#			fold=ifelse(genedata[[t]][,functypecode]==7,beta_gc[1],
#			            ifelse(genedata[[t]][,functypecode]==8,beta_gc[1]*beta_gc[2],1))
#			fe=apply(genedata[[t]][,..funcanno], 1, prod)*fold
#			pmr=fe*bmr
#			mutp[[t]]=cbind(genedata[[t]][,c(2:5)],bmr,fe,pmr)
#}
#return(mutp)
#}

#a=mutp(gene="OR4F16",funcanno = c("hic"),Afiledir = "~/SimingLab/jiezhou/driverMAPS/data",betaf0 = 2,beta_gc = c(2,2))

 #a= simulate_1funcv(gene="OR4F16",funcanno = c("hic"), Afiledir = "~/SimingLab/jiezhou/driverMAPS/data",betaf0 = 2,beta_gc = c(2,2),
          #        binary = T,para=c(0.8,0.2),Nsample=100)

#set.seed(1)
#a= power_comparebase_dev(gene="OR4F16",funcanno = c("hic"), Afiledir = "~/SimingLab/jiezhou/driverMAPS/data", binary = T,
                         betaf0 = 2,beta_gc = c(2,2), para=c(0.8,0.2),Nsample=100,Niter=2)
