#' @name oneMean
#' @aliases oneMean
#' @title CZZZ-test for one sample mean vector
#' @description Testing the equality of high dimensional mean vector to zero using the method developed in arXiv:1406.1939 [math.ST]
#'
#' @param X The \eqn{n x p}  data matrix.
#' @param m The number of Monte-Carlo samples in the test, default to be \eqn{2500}
#' @param filter A logical indicator of the filtering process, defaul to be TRUE
#' @param S Covariance matrix of \eqn{X}, if not presented it will be estimated from
#'        the input sample.
#' @param alpha The significant level of the test.
#' @param DNAME Defaul input.
#'
#' @details Implement the method developed in  arXiv:1406.1939 [math.ST] to test whether a high dimensional mean vector is zero or not, which is equivalent
#' to test \eqn{H_0: \mu=\mu_0} for some prescribed value \eqn{\mu_0} which can be subtracted from the data. The procedure utilizes bootstrap concept and derive the critical values using
#' independent Gaussian vectors whose covariance is estimated using sample covariance matrix.
#'
#' @return Value of testing statistics, p-values (the non-studentized statistic
#' and the studentized statistic respectively), alternative hypothesis, and the name of testing procedure.
#'
#' @references J. Chang, W. Zhou and W.-X. Zhou, Simulation-Based Hypothesis Testing of High Dimensional Means Under Covariance Heterogeneity (2014), arXiv:1406.1939.
#' @author Tong He
#' @export
#'
#'
#'
oneMean = function(X, m = 2500, filter = TRUE, S = NULL, alpha = 0.05, DNAME) {
  # Input Check
  checkmate::checkMatrix(X)
  checkmate::checkInt(m, lower = 1)
  checkmate::checkLogical(filter)
  # checkMatrix(S)
  checkmate::checkNumber(alpha, lower = 0, upper = 1)

  n = ncol(X)
  p = nrow(X)
  pn = 1:p

  if (filter) {
    e = 0.1
    tau1 = 0.1*(2*log(p))^(1/2-e)
    tau2 = (sqrt(2)+sqrt(2)/(2*log(p))+sqrt(2*log(1/alpha)/log(p)))*sqrt(log(p))

    # var(X,1,2) ??
    Dk = abs(rowMeans(X))/sqrt(apply(X,1,var)/n)
    pn = which(Dk >= min(tau1,tau2))
    X = X[pn,,drop = FALSE]

    p = nrow(X)
    n = ncol(X)
  }

  Tn = matrix(0, 1, 2)

  if (length(pn)>1) {
    if (!checkmate::testMatrix(S)) {
      Xm = X - rowMeans(X)
      S1 = Xm %*% t(Xm) / n
    } else {
      S1 = S[pn, pn]
    }

    Ds = diag(1/sqrt(diag(S1)))
    Sig1 = Ds %*% S1 %*% Ds

    Z3 = MASS::mvrnorm(m,rep(0,p),S1)
    Z4 = MASS::mvrnorm(m,rep(0,p),Sig1)

    Dt = abs(rowMeans(X))
    sn = diag(S1)/n
    T1 = max(sqrt(n)*Dt)
    T2 = max(Dt/sqrt(sn))

    z3.seq = apply(abs(Z3),1,max)
    z3a = quantile(z3.seq,alpha)
    z4.seq = apply(abs(Z4),1,max)
    z4a = quantile(z4.seq,alpha)
    Tn = c(T1>z3a,T2>z4a)
    names(Tn) = c("Non-Studentized", "Studentized")
    # return(Tn)
    # DNAME = deparse(substitute(X))

    names(T1) = "Non-Studentized Statistic"
    ns.p = sum(T1<z3.seq)/length(z3.seq)
    ns.res = list(statistics = T1, p.value = ns.p, alternative = "two.sided",
                  method = "One-Sample HD Non-Studentized test", data.name = DNAME)
    class(ns.res) = "htest"

    names(T2) = "Studentized Statistic"
    s.p = sum(T2<z4.seq)/length(z4.seq)
    s.res = list(statistics = T2, p.value = s.p, alternative = "two.sided",
                 method = "One-Sample HD Studentized test", data.name = DNAME)
    class(s.res) = "htest"
    res = list(NonStudent = ns.res, Student = s.res)
    return(res)

  } else {
    if (filter) {
      stop("Empty data after filtering")
    } else {
      stop("Empty data input")
    }
    return(NULL)
  }
}
