
#' @title plot phenotype, mutation and annotation for a gene across samples
#' @export
plot_mut <- function(gene_name, mut, pheno, totalnttype =96, anno_dir = ".", output_prefix = "plot", output_dir = "."){
  mycolorsori <- c("#FF3D2E","#2b8cbe","#fdae61","#41ab5d","pink","#8C2DA6","#a65628")

  afileinfo <- list(file = file.path(anno_dir, paste0("anno", totalnttype, "_nttypeXXX_annodata.txt")),
                    header = aheader,
                    coltype = acoltype,
                    totalntype = totalnttype)

  rdata <- prep_positional_data(gene_name, afileinfo, BMRreg = NULL, output_prefix = output_prefix, output_dir =  output_dir)
  fanno <-  rdata$fanno
  ri <- rdata$ri
  ri$ridx <- 1:dim(ri)[1]

  cdata <- prep_pheno_mut_data(mut, pheno)
  mut <-cdata$mut
  # add nmut info.
  nmutdt <- data.table::data.table(table(mut$SampleID))
  colnames(nmutdt) <- c("SampleID", "Nmut")
  pheno <- merge(cdata$pheno, nmutdt)
  ci <- cdata$ci

  muti <- na.omit(ci[ri[mut, on = c("chrom"= "Chromosome", "start" = "Position",  "ref" = "Ref",  "alt"= "Alt")], on = "SampleID"])
  mutmtx <- Matrix::sparseMatrix(i = muti$ridx, j = muti$cidx, dims = c(max(ri$ridx), max(ci$cidx)))
  df <-t(cbind(pheno[,2], colSums(mutmtx), colSums(mutmtx[fanno$functypecode8>0,]), colSums(mutmtx[fanno$mycons>0,]), pheno$Nmut))
  rownames(df) <- c("Phenotype", "Nonsyn", "LoF", "Cons.", "#Mut")
  df <- df[,order(df["Phenotype",])]
  par(mfrow=c(dim(df)[1],1), mar=c(0.5,5,0.5,1))
  if (length(unique(df[1,]))==2) {
    # binary phenotype
    df <- df[,order(df["#Mut",])]
    barplot(rep(1,length(df[1,])), col=c("white", "salmon")[as.factor(df[1,])], xaxt='n', yaxt = 'n', ylab= rownames(df)[1], border=NA)
  } else {
    barplot(df[1,], xaxt='n', yaxt = 'n', ylab= rownames(df)[1], border=NA)
  }
  for (i in 2:4){
    barplot(df[i,], col=mycolorsori[i], xaxt='n', yaxt = 'n', ylab= rownames(df)[i], border = mycolorsori[i])
  }
  barplot(df[5,], col=mycolorsori[5], xaxt='n', yaxt = 'n', ylim =c(0, quantile(df[5,], 0.99)[[1]]), ylab= rownames(df)[5], border=mycolorsori[5])
}
