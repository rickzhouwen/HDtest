tres = function(X) {

  # [p,n1] = size(X)
  p = nrow(X)
  n1 = ncol(X)

  # A1 = zeros(p,(n1-1)*n1)
  # B1 = zeros(p,(n1-1)*n1)
  # A2 = zeros((n1-1)*n1,p)
  # B2 = zeros((n1-1)*n1,p)
  A1 = matrix(0, p, (n1-1)*n1)
  B1 = matrix(0, p, (n1-1)*n1)
  A2 = matrix(0, (n1-1)*n1, p)
  B2 = matrix(0, (n1-1)*n1, p)

  k = 0

  for (i in 1:(n1-1)) {
    for (j in (i+1):n1) {
      k = k+1
      cat(k,'\r')
      tmp = rep(1,n1)
      tmp[c(i,j)] = 0
      mm = rowMeans(X[,which(tmp==1)])
      A1[,k] = X[,i] - mm
      A2[k,] = X[,i]
      B1[,k] = X[,j] - mm
      B2[k,] = X[,j]
    }
  }
  s1 = 0
  for (l in 1:k) {
    s1 = s1 + A2[l,] %*% B1[,l] %*% B2[l,] %*% A1[,l]
  }
  s1 = s1/(n1*(n1-1))*2
  return(s1)
}


tres12 = function(X,Y) {
  p = nrow(X)
  n1 = ncol(X)
  n2 = ncol(Y)

  mm2 = matrix(0, p ,n2)
  for (j in 1:n2) {
    tmp = rep(1, n2)
    tmp[j] = 0
    mm2[,j] = rowMeans(Y[, which(tmp==1)])
  }
  W2 = (Y - mm2) %*% t(Y)

  mm1 = matrix(0, p, n1)
  for (j in 1:n1) {
    tmp = rep(1, n1)
    tmp[j] = 0
    mm1[,j] = rowMeans(X[, which(tmp == 1)])
  }
  W1 = (X - mm1) %*% t(X)
  s12 = sum(diag(W1 %*% W2)) / (n1 * n2)
  return(s12)
}



#' @name CQ2
#' @aliases CQ2
#' @title CQ-test for two sample means
#' @description Testing the equality of two high dimensional mean vectors using the testing procedure by Chen and Qin (2010)
#'
#' @param X The n x p data matrix from the sample 1
#' @param Y The n x p data matrix from the sample 2.
#' @param DNAME Default input.
#'
#' @details Implementing testing procedure proposed by Chen and Qin (2010) to test the equality of two sample high dimensional mean vectors under the assumption of sparsity of signals.
#'
#' @return res Value of testing statistic, alternative hypothesis, and the name of testing procedure.
#'
#' @references S. Chen and Y. Qin (2010). A two-sample test for high-dimensional data with applications to gene-set testing. Ann. Statist. 38, 808-835
#' @author Tong He
#' @export
#'
#'


CQ2 = function(X,Y, DNAME) {

  p = nrow(X)
  n1 = ncol(X)
  n2 = ncol(Y)

  P1 = t(X) %*% X
  T1 = sum(P1)-sum(diag(P1))
  P2 = t(Y) %*% Y
  T2 = sum(P2)-sum(diag(P2))
  T3 = sum(sum(t(X) %*% Y))
  Tn = T1/(n1*(n1-1))+T2/(n2*(n2-1))-2*T3/(n1*n2)

  s1 = tres(X)
  s2 = tres(Y)
  s12 = tres12(X,Y)

  sn1 = 2*s1/(n1*(n1-1)) + 2*s2/(n2*(n2-1)) + 4*s12/(n1*n2)
  Qn = Tn/sqrt(sn1)
  # return(Qn)
  # DNAME = deparse(substitute(X))
  # DNAME = paste(DNAME, "and", deparse(substitute(Y)))

  names(Qn) = "Statistic"
  res = list(statistics = Qn, p.value = NA, alternative = "two.sided",
             method = "Two-Sample CQ test", data.name = DNAME)
  class(res) = "htest"
  return(res)
}

