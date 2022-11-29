#' Title
#'
#' @param b
#' @param rate.s0
#' @param mutidx
#'
#' @return
#' @export
#'
#' @examples
get_ll_s <- function(b, rate_s0, mutidx){
  rate_s <- rate_s0 * exp(b)
  rmtx <- log(1-rate_s)
  rmtx[mutidx] <- log(rate_s[mutidx])
  # log likelihood for each sample under selection
  colSums(rmtx) # faster than `colSums(log(rate.s * mut +  (1-rate.s) * (1-mut)))`
}
