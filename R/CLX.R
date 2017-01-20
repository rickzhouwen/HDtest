#' @name CLX
#' @aliases CLX
#' @title CLX-test for two sample means
#' @description Testing the equality of two high dimensional mean vectors using the testing procedure by Cai, Liu and Xia (2014).
#'
#' @param X The n x p data matrix from the sample 1
#' @param Y The n x p data matrix from the sample 2.
#' @param alpha The prescribed level of significance
#' @param DNAME Default input.
#'
#' @details Implementing testing procedure proposed by Cai, Liu, and Xia (2014) to test the equality of two sample high dimensional mean vectors under the assumption of sparsity of signals.
#'
#' @return Value of testing statistic, p-value, alternative hypothesis, and the name of testing procedure.
#'
#' @references T. Cai, W. Liu, and Y. Xia (2014). Two-sample test of high dimensional means under dependence. J. R. Statist. Soc. B. 76, 349--372
#' @author Tong He
#' @export
#'
#'
#'
#' @return
#'
CLX = function(X,Y,alpha, DNAME) {
  p = nrow(X)
  n1 = ncol(X)
  n2 = ncol(Y)

  Xbar = rowMeans(X,2)
  Ybar = rowMeans(Y,2)
  Xm = X - Xbar
  Ym = Y - Ybar

  SS = (Xm %*% t(Xm) + Ym %*% t(Ym))/(n1+n2)

  btheta = matrix(0, p, p)
  for (i in 1:p) {
    for (j in i:p) {
      xi = Xm[i,]
      xj = Xm[j,]
      yi = Ym[i,]
      yj = Ym[j,]
      btheta[i,j] = (sum((xi*xj - SS[i,j])^2) + sum((yi*yj - SS[i,j])^2))/(n1+n2)
    }
  }

  theta = btheta + t(btheta)
  theta = theta - diag(diag(theta)/2)
  Sigma = SS * as.numeric( abs(SS) >= 2*sqrt(theta*log(p)/(n1+n2)))
  eig = eigen(Sigma)
  s = diag(eig$values)
  v = eig$vectors
  Sigma = v %*% diag(pmax(diag(s), log(p)/(n1+n2))) %*% t(v)

  X = solve(Sigma, X)
  Y = solve(Sigma, Y)

  Xbar = rowMeans(X)
  Ybar = rowMeans(Y)
  Xm = X - Xbar
  Ym = Y - Ybar

  S1 = (Xm %*% t(Xm))/n1
  S2 = (Ym %*% t(Ym))/n2
  Z = Xbar - Ybar
  om = n1/(n1+n2)*diag(S1) + n2/(n1+n2)*diag(S2)
  Momega = n1*n2/(n1+n2)*max(Z^2/om)

  # DNAME = deparse(substitute(X))
  # DNAME = paste(DNAME, "and", deparse(substitute(Y)))
  Xn = as.numeric( Momega >= 2*log(p) - log(log(p)) - log(pi) - 2* log(log(1/(1-alpha))))
  p.val = exp(exp(-(Momega - 2*log(p) + log(log(p)) + log(pi))/2))
  p.val = 1- 1/p.val
  # res = c(Xn, Momega)
  names(Momega) = "Statistic"
  res = list(statistics = Momega, p.value = p.val, alternative = "two.sided",
             method = "Two-Sample CLX test", data.name = DNAME)
  class(res) = "htest"
  return(res)
}
