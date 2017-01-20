#' @name equalCovs
#' @aliases equalCovs
#' @title LC-test for equality of high dimensional covariances
#' @description Testing the equality of two high dimensional covariance matrices using the testing procedure by Li and Chen (2012).
#'
#' @param X The n x p data matrix from the sample 1
#' @param Y The n x p data matrix from the sample 2.
#' @param alpha The prescribed level of significance
#' @param DNAME Default input.
#'
#' @details Implementing testing procedure proposed by Li and Chen (2012) to test the equality of two sample high dimensional covariance matrices.
#'
#' @return Value of testing statistic, p-value, alternative hypothesis, and the name of testing procedure.
#'
#' @references J. Li and S. Chen (2012). Two sample tests for high-dimensional covariance matrices. Ann. Statist. 40, 908--940
#' @author Tong He
#' @export
#'


equalCovs <-function(X, Y, alpha, DNAME){
  # equalCovs <-function(sam1,sam2,size1,size2){
  ########################################################
  # obtain test statistic given in eqn (2.1) of the paper
  size1 = nrow(X)
  size2 = nrow(Y)
  A_mat <- X %*% t(X)
  out <- 0
  storage.mode(A_mat) <- "double"
  storage.mode(out) <- "double"
  nr <- as.integer(size1)

  # find A1
  #dyn.load("one1.dll")
  result1 <- .Fortran("code1", nr, A_mat, out=out, PACKAGE = "HDtest")
  A1 <- 2/(size1*(size1-1))*result1[[3]]

  # find A2
  #dyn.load("two2.dll")
  result2 <- .Fortran("code2", nr, A_mat, out=out, PACKAGE = "HDtest")
  A2 <- 4/(size1*(size1-1)*(size1-2))*result2[[3]]

  # find A3
  #dyn.load("three3.dll")
  result3 <- .Fortran("code3", nr, A_mat, out=out, PACKAGE = "HDtest")
  A3 <- 8/(size1*(size1-1)*(size1-2)*(size1-3))*result3[[3]]

  # obtain the statistic given by eqn (2.1) in our paper
  A_n1 <- A1-A2+A3

  # consider the sample 2

  B_mat <- Y %*% t(Y)
  out <- 0
  storage.mode(B_mat) <- "double"
  storage.mode(out) <- "double"
  nrr <- as.integer(size2)
  # find B1
  #dyn.load("one1.dll")
  result4 <- .Fortran("code1", nrr, B_mat, out=out, PACKAGE = "HDtest")
  B1 <- 2/(size2*(size2-1))*result4[[3]]

  # find B2
  #dyn.load("two2.dll")
  result5 <- .Fortran("code2", nrr, B_mat, out=out, PACKAGE = "HDtest")
  B2 <- 4/(size2*(size2-1)*(size2-2))*result5[[3]]

  # find B3
  #dyn.load("three3.dll")
  result6 <- .Fortran("code3", nrr, B_mat, out=out, PACKAGE = "HDtest")
  B3 <- 8/(size2*(size2-1)*(size2-2)*(size2-3))*result6[[3]]

  B_n2 <- B1-B2+B3

  ############################################################
  # obtain the test statistic given in eqn (2.2) of the paper
  ############################################################
  C_mat1 <- X %*% t(Y)
  C_mat2 <- Y %*% t(X)
  out <- 0
  storage.mode(C_mat1) <- "double"
  storage.mode(C_mat2) <- "double"
  storage.mode(out) <- "double"
  nrrr <- as.integer(size1)
  nl <- as.integer(size2)
  # find C1
  #dyn.load("four4.dll")
  result7 <- .Fortran("code4", nrrr, nl, C_mat1, out=out, PACKAGE = "HDtest")
  C1 <- -2/(size1*size2)*result7[[4]]

  # find C2
  #dyn.load("five5.dll")
  result8 <- .Fortran("code5", nl, nrrr, C_mat2, out=out, PACKAGE = "HDtest")
  C2 <- 4/(size1*size2*(size1-1))*result8[[4]]

  # find C3
  #dyn.load("five5.dll")
  result9 <- .Fortran("code5", nrrr, nl, C_mat1, out=out, PACKAGE = "HDtest")
  C3 <- 4/(size1*size2*(size2-1))*result9[[4]]

  # find C4
  #dyn.load("six6.dll")
  result10 <- .Fortran("code6", nrrr, nl, C_mat1, out=out, PACKAGE = "HDtest")
  C4 <- -4/(size1*size2*(size1-1)*(size2-1))*result10[[4]]

  # test statistic given by eqn (2.2) of the paper
  C_n <- C1+C2+C3+C4

  # the estimator
  T_n <- A_n1+B_n2+C_n

  # the standard deviation
  Sd_prime <- 2*(1/size1+1/size2)*((size1/(size1+size2))*A_n1+(size2/(size1+size2))*B_n2)

  test_stat <- T_n/Sd_prime
  pvalue <- 1-pnorm(test_stat)
  decision <- as.numeric(pvalue < alpha)
  # test <- matrix(0, 3, 1)
  # test[,1] <- c(decision, pvalue, test_stat)
  # colnames(test) = "LC"
  # rownames(test) = c('decision', 'p-value', 'test statistic')
  # return(test)
  # DNAME = deparse(substitute(X))
  # DNAME = paste(DNAME, "and", deparse(substitute(Y)))

  names(test_stat) = "Statistic"
  res = list(statistics = test_stat, p.value = pvalue, alternative = "two.sided",
             method = "Two-Sample LC test", data.name = DNAME)
  class(res) = "htest"

  return(res)
}


