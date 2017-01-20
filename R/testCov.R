#' @name testCov
#' @aliases testCov
#' @title Testing the equality of two sample covariance matrices in high dimension.
#' @description Testing the equality of two sample covariance matrices in high dimension using different methods.
#'
#' @param X the n x p training data, could be a \code{matrix} or a \code{data.frame} object.
#' @param Y the n x p training data matrix, could be a \code{matrix} or a \code{data.frame} object.
#' @param method a string incidating the method for the test. The current available
#'        methods are \code{ALL}, \code{HD}, \code{LC}, \code{CLX}, \code{Scott}.
#' @param J the number of repetition in the test
#' @param alpha the significant level of the test.
#' @param n.core the number of cores to be used in parallel when \code{HD}
#'        is called.
#'
#'
#' @return
#'
#' For any single method, the function returns an \code{htest} object.
#'
#' For method \code{ALL}: A list of four \code{htest} objects.
#'
#' HD refers to "Chang, J., Zhou, W., Zhou, W.-X., and Wang, L. (2016). Comparing large covariance
#' matrices under weak conditions on the dependence structure and its application to gene
#'  clustering. Biometrics. To appear"#'
#'
#' CLX refers to "Cai, T. T., Liu, W., and Xia, Y. (2013).
#' Two-sample covariance matrix testing and support recovery in high-dimensional
#' and sparse settings. Journal of the American Statistical Association 108, 265-277."
#'
#' Sc refers to "Schott, J. R. (2007). A test for the equality of covariance
#' matrices when the dimension is large relative to the sample size.
#' Computational Statistics and Data Analysis 51, 6535-6542."
#'
#' @examples
#' data(GO54)
#' testCov(GO54$X, GO54$Y, method = "ALL", J = 100)
#' data(GO26)
#' testCov(GO26$X, GO26$Y, method = "ALL", J = 100)
#'
#' @author Tong He
#' @export
#'
#'
testCov = function(X, Y, method = "ALL", J = 2500, alpha = 0.05, n.core = 1) {
  DNAME = deparse(substitute(X))
  DNAME = paste(DNAME, "and", deparse(substitute(Y)))
  X = data.matrix(X)
  Y = data.matrix(Y)
  if (method == "HD") {
    res = testCovHD(X = X, Y = Y, J = J, alpha = alpha, DNAME = DNAME,
                    n.core = n.core)
    res = res[[1]]
  } else if (method == "LC") {
    res = equalCovs(X = X, Y = Y, DNAME = DNAME)
  } else if (method == "CLX") {
    res = testCovHD(X = X, Y = Y, J = 1, alpha = alpha, DNAME = DNAME,
                    n.core = 1)
    res = res[[2]]
  } else if (method == "Scott") {
    message("The method \'Scott\' does not support skewed data.")
    res = testCovHD(X = X, Y = Y, J = 1, alpha = alpha, DNAME = DNAME,
                    n.core = 1)
    res = res[[3]]
  } else if (method == "ALL") {
    res = testCovHD(X = X, Y = Y, J = J, alpha = alpha, DNAME = DNAME,
                    n.core = n.core)
    # nms = colnames(res)
    res2 = equalCovs(X = X, Y = Y, alpha = alpha, DNAME = DNAME)
    res = c(res, list(res2))
    # colnames(res) = c(nms, 'LC')
    message("The method \'Scott\' does not support skewed data.")
  } else {
    msg = paste0("method \'", method, "\' not available.")
    stop(msg)
  }
  return(res)
}


testCovHD = function(X, Y, J = 2500, alpha = 0.05, DNAME, n.core = 1) {
  checkmate::checkMatrix(X)
  checkmate::checkMatrix(Y)

  p = ncol(X)
  n1 = nrow(X)
  n2 = nrow(Y)
  if (ncol(Y) != p) {
    stop('Different dimensions of X and Y.')
  }

  tmp = rep(c(rep(1, n1)/n1, rep(1, n2)/n2), J)
  scalev = matrix(tmp, ncol = 1)

  qalpha = -log(8*pi) - 2*log(log(1/(1-alpha)))
  cri = 4*log (p)-log (log (p)) + qalpha


  Sx = cov(X)*(n1-1)/n1
  Sy = cov(Y)*(n2-1)/n2

  xa = t(t(X) - colMeans(X))
  ya = t(t(Y) - colMeans(Y))

  vx = t(xa^2)%*%(xa^2)/n1 - 2/n1 * (t(xa)%*% xa) * Sx + Sx^2
  vy = t(ya^2)%*%(ya^2)/n2 - 2/n2 * (t(ya)%*% ya) * Sy + Sy^2

  numo = abs(Sx-Sy)
  deno = sqrt(vx/n1 + vy/n2)
  Tnm = max(numo/deno)

  xat = t(xa)/n1
  yat = t(ya)/n2
  # g = rnorm((n1+n2)*J)*scalev

  # ts = matrix(0,J,1)
  scalev = c(rep(1, n1)/n1, rep(1, n2)/n2)
  if (n.core > 1) {
    # for (j in 1:J) {
    cl = makeCluster(n.core)
    registerDoParallel(cl)
    ts = foreach(j = 1:J, .combine = rbind) %dopar% {
      # ind1 = ((j-1)*(n1+n2)+1):((j-1)*(n1+n2)+n1)
      # ind2 = ((j-1)*(n1+n2)+n1+1):((j-1)*(n1+n2)+n2+n1)

      g = rnorm(n1+n2)*scalev
      atmp = sum(g[1:n1])
      btmp = sum(g[(n1+1):(n1+n2)])

      ts1 = (t(xa*g[1:n1]) - xat*atmp)%*% xa
      ts2 = (t(ya*g[(n1+1):(n1+n2)]) - yat*btmp)%*% ya

      # atmp = sum(g[ind1])
      # btmp = sum(g[ind2])

      # ts1 = (t(xa*g[ind1]) - xat*atmp)%*% xa
      # ts2 = (t(ya*g[ind2]) - yat*btmp)%*% ya
      #ts[j] = max(abs(ts1-ts2)/deno)
      max(abs(ts1-ts2)/deno)
    }
    stopCluster(cl)
  } else {
    ts = matrix(0,J,1)
    for (j in 1:J) {
      g = rnorm(n1+n2)*scalev
      atmp = sum(g[1:n1])
      btmp = sum(g[(n1+1):(n1+n2)])

      ts1 = (t(xa*g[1:n1]) - xat*atmp)%*% xa
      ts2 = (t(ya*g[(n1+1):(n1+n2)]) - yat*btmp)%*% ya

      # ind1 = ((j-1)*(n1+n2)+1):((j-1)*(n1+n2)+n1)
      # ind2 = ((j-1)*(n1+n2)+n1+1):((j-1)*(n1+n2)+n2+n1)
      #
      # atmp = sum(g[ind1])
      # btmp = sum(g[ind2])
      #
      # ts1 = (t(xa*g[ind1]) - xat*atmp)%*% xa
      # ts2 = (t(ya*g[ind2]) - yat*btmp)%*% ya
      ts[j] = max(abs(ts1-ts2)/deno)
    }
  }

  ZCZt = Tnm > quantile(ts, 1-alpha)
  CLX = max((Sx-Sy)^2/(vx/n1+vy/n2))
  # CLX = CLX-(4*log(p)-log(log(p)))

  Sxx = Sx*n1/(n1-1)
  Syy = Sy*n2/(n2-1)

  SsS = (Sxx*n1 + Syy*n2)/(n1+n2)
  eta1 = ((n1-1)+2)*((n1-1)-1)
  eta2 = ((n2-1)+2)*((n2-1)-1)
  d1 = (1-(n1-1-2)/eta1)*sum(diag(Sxx%*%Sxx))
  d2 = (1-(n2-1-2)/eta2)*sum(diag(Syy%*%Syy))
  d3 = 2*sum(diag(Sxx %*% Syy))
  d4 = (n1-1)/eta1*sum(diag(Sxx))^2
  d5 = (n2-1)/eta2*sum(diag(Syy))^2
  th = 4*(((n1+n2-2)/((n1-1)*(n2-1)))^2)*
    ((n1+n2-2)^2/((n1+n2)*(n1+n2-2-1))*
       (sum(diag(SsS %*% SsS))-(sum(diag(SsS)))^2/(n1+n2-2)))^2
  Sc = (d1+d2-d3-d4-d5)/sqrt(th)

  # Res = matrix(0,3,3)
  # Res[1, ] = c(ZCZt, CLX>cri, abs(Sc)>qnorm(1-alpha/2))
  # Res[3, ] = c(Tnm, CLX, abs(Sc))
  # CLX = CLX-(4*log(p)-log(log(p)))
  # Res[2, ] = c(sum(ts>=Tnm)/J,
  #              1 - exp(-exp(-CLX/2)/(sqrt(8*pi))),
  #              (1-pnorm(abs(Sc)))*2)
  #
  # rownames(Res) = c('decision', 'p-value', 'test statistic')
  # colnames(Res) = c('HDtest', 'CLX', 'Scott')
  # return(Res)

  # DNAME = deparse(substitute(X))
  # DNAME = paste(DNAME, "and", deparse(substitute(Y)))

  names(Tnm) = "Statistic"
  hd.res = list(statistics = Tnm, p.value = sum(ts>=Tnm)/J, alternative = "two.sided",
                method = "Two-Sample HD test", data.name = DNAME)
  class(hd.res) = "htest"

  names(CLX) = "Statistic"
  CLX.rev = CLX-(4*log(p)-log(log(p)))
  clx.p = 1 - exp(-exp(-CLX.rev/2)/(sqrt(8*pi)))
  clx.res = list(statistics = CLX, p.value = clx.p, alternative = "two.sided",
                 method = "Two-Sample CLX test", data.name = DNAME)
  class(clx.res) = "htest"

  names(Sc) = "Statistic"
  sc.p = (1-pnorm(abs(Sc)))*2
  sc.res = list(statistics = abs(Sc), p.value = sc.p, alternative = "two.sided",
                method = "Two-Sample Scott test", data.name = DNAME)
  class(sc.res) = "htest"

  res = list(HD = hd.res, CLX = clx.res, Scott = sc.res)
  return(res)
}


