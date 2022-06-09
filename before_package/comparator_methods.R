mlr <- function(mutdatalist,edata){
  mutdata <- do.call(rbind, mutdatalist)
  dstatus <- colSums(mutdata)
  dstatus[which(dstatus>0)] <- 1
  mlrfit <- lm(edata~dstatus)
  mlrres <- summary(mlrfit) 
  return(mlrres)
}

genebinom <- function(mutdatalist,edata){
  # binom on gene level
  gmut <- do.call(rbind,mutdatalist)
  mu.c <- length(edata[edata==1])
  mu.n <- length(edata[edata==0])
  gmutc <- sum(gmut[, 1:mu.c])
  gmutn <- sum(gmut[, (mu.c+1): dim(gmut)[2]])
  resp <- binom.test(gmutc, gmutn+gmutc, p=mu.c/(mu.c+mu.n))$p.value
  return(resp)
}

genelr <- function(mutdatalist,edata){
  # logistic regression on gene level
  mutdata <- do.call(rbind, mutdatalist)
  dstatus <- colSums(mutdata)
  dstatus[which(dstatus>0)] <- 1
  glrfit <- glm(edata~dstatus,family = binomial(link = "logit"))
  glrres <- summary(glrfit)
  return(glrres)
}

genefisher <- function(mutdatalist,edata){
  # fisher's test on gene level
  gmut <- do.call(rbind,mutdatalist)
  mu.c <- length(edata[edata==1])
  mu.n <- length(edata[edata==0])
  gmutc <- sum(gmut[, 1:mu.c])
  gmutn <- sum(gmut[, (mu.c+1): dim(gmut)[2]])
  resp <- fisher.test(matrix(c(gmutc, gmutn,mu.c, mu.n), ncol=2))$p.value
  return(resp)
}
