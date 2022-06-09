#' Estimate the effect of sample level variables on background mutation rate
#' @param canno A matix with column Phenotype indicating the number of the silence mutations for every sample
#' @param ganno A matrix with column baseline giving position level bmr.
#' @param l A vector indicating the lengths of the genes involved
#' @param beta A vector giving the coefficients of position level variables
#' @param X A data frame giving position level variables
#' @param S A data frame giving sample level variables
#' @return A list giving all the information in the model
bmr_omega=function(canno,ganno,l,beta,X,S){
# number of mutations on sample level.
Y=canno$Phenotype
baseline=ganno$baseline
# sample size
n=length(Y)
Z=exp(X%*%beta)
U=c()
ll=cumsum(l)
U[1]=sum(Z[1:ll[1]])
for (i in 2:m) {
U[i]=sum(Z[(ll[i-1]+1):ll[i]])
}
offset=log(sum(U*lambda))
preg=glm(Y~S,offset = offset,family = "poisson")
logmu=pref$coefficients[1,]
omega=preg$coefficients[2,]
muij=vector("list",n)

for (i in 1:n) {
  muij[[i]]=c()
  for (j in 1:ll[1]) {
    muij[[i]]=c(muij[[i]],exp(logmu[1])*exp(X[j,]%*%beta)*exp(S[i,]%*%omega))
  }
for (g in 2:m) {
  for (j in ll[g-1]:ll[g]) {
muij[[i]]=c(muij[[i]],exp(logmu[1])*exp(X[j,]%*%beta)*exp(S[i,]%*%omega))
  }
}
}
return(list(logmu,omega,muij))
}

