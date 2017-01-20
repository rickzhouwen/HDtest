#' @name segment
#' @aliases segment
#' @title segment
#' @description segment
#'
#' @param Y \eqn{p x n} data matrix, with \eqn{p} time series of length \eqn{n}
#' @param mean_y mean vector
#' @param k lag
#' @param n each time series has length \eqn{n}
#' @param p \eqn{p} time series
#'
#' @return b An optimal bandwidth.
#'
#' @references J. Chang, Q. Yao, and W. Zhou (2016) Testing for high-dimensional white noise using maximum cross correlations. \emph{Biometrika}, to appear.
#' @author Meng Cao
#' @export
#'
segment = function(Y, mean_y, k, n, p){
  b = rep(0, p^2)
  for (t in (1:(n-k))){
    s = 0
    C = Y[,(t+k)] - mean_y
    D = Y[, t] - mean_y
    for (i in (1:p)){
      for (j in (1:p)){
        s = s+1
        b[s] = b[s] + C[i]*D[j]
      }
    }
  }
  b = b/n
  return (b)
}
