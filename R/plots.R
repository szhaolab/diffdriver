#' @title plot mutations
#' @export
plot_mut <- function(mutmtx, canno, bmrmtx, ganno){
  mycolorsori <- c("#FF3D2E","#2b8cbe","#fdae61","#41ab5d","pink","#8C2DA6","#a65628")

  df <- rbind(canno$Phenotype,colSums(mutmtx), colSums(mutmtx[ganno$functypecode8>0,]), colSums(mutmtx[ganno$mycons>0,]),colSums(exp(bmrmtx)))
  rownames(df) <- c("Phenotype", "Nonsyn", "LoF", "Cons.", "BMR")
  df <- df[,order(df["Phenotype",])]
  par(mfrow=c(dim(df)[1],1), mar=c(0.5,5,0.5,1))
  if (unique(df[1,])==2) {
    # binary phenotype
    df <- df[,order(df["BMR",])]
    barplot(rep(1,length(df[1,])), col=c("white", "salmon")[as.factor(df[1,])], xaxt='n', yaxt = 'n', ylab= rownames(df)[1], border=NA)
  } else {
    barplot(df[1,], xaxt='n', yaxt = 'n', ylab= rownames(df)[1], border=NA)
  }
  for (i in 2:4){
    barplot(df[i,], col=mycolorsori[i], xaxt='n', yaxt = 'n', ylab= rownames(df)[i], border = mycolorsori[i])
  }
  barplot(df[5,], col=mycolorsori[5], xaxt='n', yaxt = 'n', ylim =c(0, quantile(df[5,], 0.99)[[1]]), ylab= rownames(df)[5], border=mycolorsori[5])
}
