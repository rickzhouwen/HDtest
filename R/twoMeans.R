#' @name twoMeans
#' @aliases twoMeans
#' @title CZZZ-test for two sample mean vectors
#' @description Testing the equality of two sample high dimensional mean vectors using the method developed in arXiv:1406.1939 [math.ST]
#'
#'
#' @param X The n x p training data matrix.
#' @param Y The n x p training data matrix.
#' @param m The number of repetition in the test, default to be 2500
#' @param filter A logical indicator of the filtering process, default to be TRUE
#' @param SX The covariance matrix of X, if not presented it will be estimated from
#'        the input sample.
#' @param SY The covariance matrix of T, if not presented it will be estimated from
#'        the input sample.
#' @param alpha The significant level of the test.
#' @param DNAME Defaulf input.
#'
#' @details Implement the method developed in  arXiv:1406.1939 [math.ST] to test whether a high dimensional mean vector is zero or not, which is equivalent
#' to test \eqn{H_0: \mu_1=\mu_2}. The procedure utilizes bootstrap concept and derive the critical values using
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
#'
twoMeans = function(X, Y, m = 2500, filter = TRUE, SX = NULL, SY = NULL, alpha = 0.05,
                    DNAME) {
  checkmate::checkMatrix(X)
  checkmate::checkMatrix(Y)
  checkmate::checkInt(m, lower = 1)
  checkmate::checkLogical(filter)
  # checkMatrix(S)
  checkmate::checkNumber(alpha, lower = 0, upper = 1)

  p = nrow(X)
  pn = 1:p
  n1 = ncol(X)
  n2 = ncol(Y)
  if (nrow(Y) != p) {
    stop('Different dimensions of X and Y.')
  }

  if (filter) {
    e = 0.1
    tau1 = 0.1*(2*log(p))^(1/2-e)
    tau2 = (sqrt(2)+sqrt(2)/(2*log(p))+sqrt(2*log(1/alpha)/log(p)))*sqrt(log(p))

    Dk = abs(rowMeans(X) - rowMeans(Y))/sqrt(apply(X,1,var)/n1+apply(Y,1,var)/n2)
    pn = which(Dk >= min(tau1,tau2))
    X = X[pn,,drop = FALSE]
    Y = Y[pn,,drop = FALSE]

    p = nrow(X)
    n1 = ncol(X)
    n2 = ncol(Y)
  }

  Tn = matrix(0, 1, 2)

  if (length(pn)>1) {
    if (!checkmate::testMatrix(SX) || !checkmate::testMatrix(SY)) {
      Xm = X - rowMeans(X)
      S1 = Xm %*% t(Xm) / n1

      Ym = Y - rowMeans(Y)
      S2 = Ym %*% t(Ym) / n2
    } else {
      S1 = SX[pn, pn]
      S2 = SY[pn, pn]
    }

    Sig1 = S1+n1/n2*S2
    Ds = diag(1/sqrt(diag(Sig1)))
    Sig2 = Ds %*% Sig1 %*% Ds

    Z3 = MASS::mvrnorm(m,rep(0,p),Sig1)
    Z4 = MASS::mvrnorm(m,rep(0,p),Sig2)

    #     z3a = quantile(apply(abs(Z3),1,max),alpha)
    #     z4a = quantile(apply(abs(Z4),1,max),alpha)

    Dt = abs(rowMeans(X) - rowMeans(Y))
    sn = diag(Sig1)/n1
    T1 = max(sqrt(n1)*Dt)
    T2 = max(Dt/sqrt(sn))

    # z3a = quantile(apply(abs(Z3),1,max),alpha)
    # z4a = quantile(apply(abs(Z4),1,max),alpha)
    z3.seq = apply(abs(Z3),1,max)
    z3a = quantile(z3.seq,alpha)
    z4.seq = apply(abs(Z4),1,max)
    z4a = quantile(z4.seq,alpha)
    Tn = c(T1>z3a,T2>z4a)
    names(Tn) = c("Non-Studentized", "Studentized")
    # return(Tn)
    # DNAME = deparse(substitute(X))
    # DNAME = paste(DNAME, "and", deparse(substitute(Y)))

    names(T1) = "Non-Studentized Statistic"
    ns.p = sum(T1<z3.seq)/length(z3.seq)
    ns.res = list(statistics = T1, p.value = ns.p, alternative = "two.sided",
                  method = "Two-Sample HD Non-Studentized test", data.name = DNAME)
    class(ns.res) = "htest"

    names(T2) = "Studentized Statistic"
    s.p = sum(T2<z4.seq)/length(z4.seq)
    s.res = list(statistics = T2, p.value = s.p, alternative = "two.sided",
                 method = "Two-Sample HD Studentized test", data.name = DNAME)
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
