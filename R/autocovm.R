#' @name autocovm
#' @aliases autocovm
#' @title Lag-k autocovariance matrix for multivariate time series
#' @description This is an auxillary function to compute the autocovariance matrix for multivariate time series at lag k.
#' @param Y A multivariate time series.
#' @param k The lag k.
#'
#' @details  Compute the autocovariance matrix of mutlivariate time series \eqn{Y} at lag \eqn{k}.
#' @return sm The autocovariace matrix at lag k.
#' @import MASS
#' @author Meng Cao, Wen Zhou
#'
#' @export
autocovm = function(Y, k){
  n = dim(Y)[1]
  p = dim(Y)[2]
  Y = t(Y)
  sm = rep(0, p)
  a1 = t(t(Y[, 1:(n-k)]))
  a2 = t(t(Y[, (1+k):n]))
  a1x = sweep(a1, MARGIN = 1, (apply(X = a1, MARGIN = 1, FUN = mean)), FUN = "-")
  a2x = sweep(a2, MARGIN = 1, (apply(X = a2, MARGIN = 1, FUN = mean)), FUN = "-")
  for (t in c(1:(n-k))){
    sm = sm + a1x[, t] %*% t(a2x[, t])
  }
  sm = sm/(n-k-1)
  return(sm)
}
