#' @name aeS
#' @aliases aeS
#' @title Quantile of the absolute values of Gaussian vectors with long run covariance
#' @description This is an auxillary function to compute alpha-level quantile of absolute values of Gaussian vectors whose covariance matrices are specified by W.
#'
#' @param ft The input multivariate time series
#' @param Sn Long run covariance matrix
#' @param W The inverse of covariance matrix at lag 0
#' @param M Number of Gaussian vectors sampled for computing empirical quantile
#' @param alpha level of significance (default 0.05)
#' @details For input multivariate time series \eqn{\varepsilon_t}, derive
#' \deqn{f_t=\{vec(\varepsilon_{t+1}\varepsilon_t^T),\ldots,vec(\varepsilon_{t+K}\varepsilon_t^T)\}^T}. Long run covariance matrix
#' \eqn{S_n} is estimated separately using the method described in Section 2.3 in the reference below and inverse
#' covariance matrix at lag \eqn{0} is estimated using \eqn{\varepsilon_t}. \eqn{M} independent Gaussian vectors with desired long run covariance are sampled to
#' compute the \eqn{\alpha}-level empirical quantiles for their absolute values.
#'
#' @return alpha-level quantiles for testing whether the input multivariate time series is a white noise or not, alpha is 0.05
#'
#' @references J. Chang, Q. Yao, and W. Zhou (2016) Testing for high-dimensional white noise using maximum cross correlations. \emph{Biometrika}, to appear.
#' @author Meng Cao, Wen Zhou
#' @export
#'

aeS = function(ft, Sn, W, M, alpha){
  p = dim(ft)[1]
  n = dim(ft)[2]
  xi = t(mvrnorm(n = M, mu = rep(0, n), Sn))/sqrt(n)
  G1 = ft %*% xi
  G = W*G1
  cv = quantile(apply(abs(G),2, max), 1-alpha)
  return(cv)
}

