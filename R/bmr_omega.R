#' Estimate the effect of sample level variables on background mutation rate
#' @param Y is a vector indicating the number of the silence mutations for every subject
#' @param lambda A vector giving  gene effects
#' @param l A vector indicating the lengths of the genes involved
#' @param beta A vector giving the coefficients of position level variables
#' @param X A data frame giving position level variables
#' @param S A data frame giving sample level variables
#' @return A list giving all the information in the model
bmr_omega=function(Y,lambda,l,beta,X,S){
m=length(lambda)
Z=exp(X%*%beta)
U=c()
ll=cumsum(l)
U[1]=sum(Z[1:ll[1]])
for (i in 2:m) {
U[i]=sum(Z[(ll[i-1]+1):ll[i]])
}
offset=log(sum(U*lambda))
preg=glm(Y~S,offset = offset)
omega=preg$coefficients[2,]
return(omega)
}
