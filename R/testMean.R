#' @name testMean
#' @aliases testMean
#' @title Testing the equality of two sample mean vectors in high dimension.
#' @description Testing the equality of two sample mean vectors in high dimension using different methods.
#'
#'
#' @param X the n x p training data matrix, could be a \code{matrix} or a \code{data.frame} object.
#' @param Y the n x p training data matrix, if presented the method will perform
#'          a two-sample test of mean, one-sample test otherwise.
#'          Could be a \code{matrix} or a \code{data.frame} object.
#' @param method a string incidating the method for the test. The current available
#'        methods are \code{ALL}, \code{HD}, \code{CQ}, \code{CLX}.
#' @param m the number of repetition in the test
#' @param filter a logical indicator of the filtering process
#' @param alpha the significant level of the test.
#' @param SX covariance matrix of X, if not presented it will be estimated from
#'        the input sample.
#' @param SY covariance matrix of T, if not presented it will be estimated from
#'        the input sample.
#'
#' @return
#' For method \code{HD}, the function returns two \code{htest} objects for
#' non-studentized and studentized test respectively.
#'
#' For method \code{CLX} and \code{CQ}, the function returns an \code{htest} object.
#'
#' For method \code{ALL}: A list of four \code{htest} objects.
#'
#' HD refers to arXiv:1406.1939 [math.ST]
#'
#' @examples
#' data(GO54)
#' testMean(GO54$X, m = 100, method = "HD")
#' testMean(GO54$X, GO54$Y, m = 100, method = "ALL")
#'
#' @author Tong He
#' @export
#'
testMean = function(X, Y = NULL, method = "HD", m = 2500, filter = TRUE,
                    alpha = 0.05, SX = NULL, SY = NULL) {
  DNAME = deparse(substitute(X))
  X = data.matrix(X)
  if (!is.null(Y)) {
    DNAME = paste(DNAME, "and", deparse(substitute(Y)))
    Y = data.matrix(Y)
  }
  if (is.null(Y)) {
    # DNAME = deparse(substitute(X))
    if (method != "HD") {
      stop("One sample test is only supported by method \'HD\'.")
    }
    res = oneMean(X = t(X), m = m, filter = filter, S = SX, alpha = alpha,
                  DNAME = DNAME)
  } else {
    # DNAME = deparse(substitute(X))
    # DNAME = paste(DNAME, "and", deparse(substitute(Y)))
    if (method == "HD") {
      if (!is.null(SY)) {
        stop("Y is not presented, thus SY should not be presented as well.")
      }
      res = twoMeans(X = t(X), Y = t(Y), m = m, filter = filter,
                     SX = SX, SY = SY, alpha = alpha, DNAME = DNAME)
    } else if (method == "CLX") {
      res = CLX(X = t(X), Y = t(Y), alpha = alpha, DNAME = DNAME)
    } else if (method == "CQ") {
      res = CQ2(X = t(X), Y = t(Y), DNAME = DNAME)
    } else if (method == "ALL") {
      res1 = twoMeans(X = t(X), Y = t(Y), m = m, filter = filter,
                      SX = SX, SY = SY, alpha = alpha, DNAME = DNAME)
      res2 = CLX(X = t(X), Y = t(Y), alpha = alpha, DNAME = DNAME)
      res3 = CQ2(X = t(X), Y = t(Y), DNAME = DNAME)
      res = list(HD = res1, CLX = res2, CQ = res3)
    } else {
      msg = paste0("method \'", method, "\' not available.")
      stop(msg)
    }
  }
  return(res)
}
