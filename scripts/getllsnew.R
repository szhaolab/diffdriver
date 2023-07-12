#' Title
#'
#' @param rate.s0
#' @return
#' @export
#'
#' @examples
get_ll_s <- function(rate_s0,b){
  # log likelihood for each sample under selection
if (any(rate_s0<=0)){stop("rate_s0 should be positive!")}
rate_s0=rate_s0*exp(b)
  colSums(log(rate_s0)) # faster than `colSums(log(rate.s * mut +  (1-rate.s) * (1-mut)))`
}
