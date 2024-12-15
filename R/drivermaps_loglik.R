# library("reshape2")
# global variables: Totalnttype

loglik1_sub <- function(vbetasub, data, y){
  lastidx <- length(vbetasub)
  y * (as.matrix(data) %*% vbetasub[1:(lastidx-1)] + vbetasub[lastidx])
}

loglik1sum <- function(vbeta, matrixlist){
  ll1 <- 0
  for (j in seq(1:Totalnttype)){
    vbetasub <- convertbeta(j, vbeta)
    anno <- matrixlist[[j]][[1]]
    y <- matrixlist[[j]][[2]]
    ddmyn  <- anno[which(y>0)]
    if (dim(ddmyn)[1] == 0) next
    yn<- y[which(y>0)]
    ll1 <- ll1 + sum(loglik1_sub(vbetasub, ddmyn, yn))
  }
  ll1
}

loglik2 <- function(alpha, y_g_s, mu_g_s){
  # each gene has a prior
  (alpha + y_g_s)*log(alpha + mu_g_s)- lgamma(alpha + y_g_s)
}

loglik3 <- function(alpha, y_g, y_g_s){
  # y_g a vector of numbers of mutations in genes
  # y_g_s a vector of numbers of silent mutations in genes (prior data)
  lgamma(alpha + y_g + y_g_s)
}

loglik4 <- function(alpha, y_g, mu_g, y_g_s, mu_g_s){
  # mu_g is a function of beta for each gene
  # mu_g_s is a function of BM related beta for each gene (expected silent mutation number)
  (alpha + y_g + y_g_s) * log(alpha + mu_g + mu_g_s)
}

gene_mu_sub <- function(vbetasub, data, gene){
  # generate mu_g
  lastidx <- length(vbetasub)
  mu <- exp(as.matrix(data) %*% vbetasub[1:(lastidx-1)] + vbetasub[lastidx])
  mu <- data.table(mu)
  mu[,agg_var:= gene]
  aggregate_bysum(mu)
}

gene_mu <- function(vbeta,matrixlist){
  mu_g <- data.table(agg_var = character(0), V1 = numeric(0))
  for (j in seq(1:Totalnttype)){
    vbetasub <- convertbeta(j, vbeta)
    ddm <- matrixlist[[j]][[1]]
    gene <- matrixlist[[j]][[3]]
    mu_g <- rbind(mu_g, gene_mu_sub(vbetasub, ddm, gene))
  }
  mu_g <- aggregate_bysum(mu_g)
  mu_g
}

gene_y <- function(matrixlist){
  # generate y_g
  y_g <- data.table(agg_var = character(0), y = numeric(0))
  for (j in seq(1:Totalnttype)){
    y <- matrixlist[[j]][[2]]
    gene <- matrixlist[[j]][[3]]
    y<- data.table(y)
    y$agg_var <- gene
    y_g <- rbind(y_g, aggregate_bysum(y))
  }
  aggregate_bysum(y_g)
}


# pars in loglikfn: 1-Totalnttype are beta_0t, last par is alpha,
# the second to last is beta_f0
loglikfn <- function(pars, matrixlist, y_g_s_in, mu_g_s_in){
  vbeta <-pars[1: length(pars)-1]
  alpha <- pars[length(pars)]
  y_g <- gene_y(matrixlist)
  mu_g <- gene_mu(vbeta,matrixlist)
  y_g_s <- y_g_s_in[y_g][,1:2, with=F]
  y_g_s[is.na(y_g_s)] <- 0
  mu_g_s <- mu_g_s_in[y_g][,1:2, with=F]
  mu_g_s[is.na(mu_g_s)] <- 0
  ll <- loglik1sum(vbeta, matrixlist) + sum(loglik2(alpha, y_g_s$y, mu_g_s$V1)) +
    sum(loglik3(alpha, y_g$y, y_g_s$y)) - sum(loglik4(alpha, y_g$y, mu_g$V1, y_g_s$y, mu_g_s$V1))
  return(ll)
}

aggregate_bysum <- function(dat, xs="agg_var") {
  # Convert to data.table.
  # dat <- data.table(dat)
  # Append the vector of group names as an extra column.
  # dat$agg_var <- xs
  # Melt the data.table so all values are in one column called "value".
  dat <- melt(dat, id.vars = xs)
  # Cast the data.table back into the original shape, and take the sum.
  dat <- dcast.data.table(
    dat, paste(xs, "~ variable", sep=""), value.var = "value",
    fun.aggregate = sum, na.rm = TRUE
  )
  #rownames(dat) <- dat$agg_var
  # Delete the extra column.
  #dat[ , agg_var := NULL]
  dat
}

##' optifix. Optimise with fixed parameters
##'
##' its like optim, but with fixed parameters.
##'
##' specify a second argument 'fixed', a vector of TRUE/FALSE values.
##' If TRUE, the corresponding parameter in fn() is fixed. Otherwise its
##' variable and optimised over.
##'
##' The return thing is the return thing from optim() but with a couple of extra
##' bits - a vector of all the parameters and a vector copy of the 'fixed' argument.
##'
##' Written by Barry Rowlingson <b.rowlingson@lancaster.ac.uk> October 2011
##'
##' This file released under a CC By-SA license:
##' http://creativecommons.org/licenses/by-sa/3.0/
##'
##' and must retain the text: "Originally written by Barry Rowlingson" in comments.
##'
optifix <- function(par, fixed, fn, gr = NULL, ...,
                    method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
                    lower = -Inf, upper = Inf,
                    control = list(), hessian = FALSE){
  force(fn)
  force(fixed)
  .npar=length(par)
  .fixValues = par[fixed]

  .parStart = par[!fixed]

  .fn <- function(par,...){
    .par = rep(NA,sum(!fixed))
    .par[!fixed] = par
    .par[fixed] = .fixValues
    fn(.par,...)
  }

  if(!is.null(gr)){
    .gr <- function(par,...){
      .gpar = rep(NA,sum(!fixed))
      .gpar[!fixed] = par
      .gpar[fixed] = .fixValues
      gr(.gpar,...)[!fixed]
    }
  }else{
    .gr <- NULL
  }

  .opt = optim(.parStart,.fn,.gr,...,method=method,lower=lower,control=control,hessian=hessian)

  .opt$fullpars = rep(NA,sum(!fixed))
  .opt$fullpars[fixed]=.fixValues
  .opt$fullpars[!fixed]=.opt$par
  .opt$fixed = fixed
  return(.opt)
}



