#' Test_our_new
#' @description Our new method to test white noise
#' @param Y is input multivaraite time series data
#' @param k_max is a parameter (for example default 10)
#' @param kk is a vector of parameters, could be a scalar as well (kk = c(2:10))
#' @param M is a parameter, could be 1000, 2000 for example
#' @param k0 is parameter in time series PCA for transformation (default 10)
#' @param delta is 2nd parameter in time series PCA for transformation (default 1.5)
#' @param type 1: wntest, 2: test_LM, 3: test_pre, 4: test_TB
#' type = 1 need X, k_max, ,kk, M, bw, type = 2 need Y, k, type = 3: need Y, k_max, kk, type = 4: need Y
#' @param alpha level of significance
#' @param opt = 1, perform transformation, else do not perform transformation
#' @author Meng Cao
#' @return res white noise or not
#' @examples
#' library(expm)
#' p = 15
#' n = 300
#' S1 = diag(1, p, p)
#' for(ii in c(1:p)){
#' for(jj in c(1:p)){
#' S1[ii, jj] = 0.995^(abs(ii-jj))
#' }
#' }
#' S11 = sqrtm(S1)
#' X = S11 %*% matrix(rt(n*p, df = 8), ncol = n)
#' k_max = 10
#' kk = seq(2, k_max, 2)
#' M = 2000
#' k0 = 10
#' delta = 1.5
#' alpha = 0.05
#' wntest(X, M, k_max, kk, type = 1, opt = 0)
#'
#' @export
wntest = function(Y, M, k_max = 10, kk, type = 1, alpha = 0.05, k0 = 10, delta = 1.5, opt = 1){
  if (type == 1){
    X = Y

    if(opt == 1){
      X = t(sgTs_thre(t(X), k0 = k0, delta = delta))
    }
    bw = opbw(X)
    p = dim(X)[1]
    n = dim(X)[2]
    R0 = diag(1/sqrt(diag(X%*%t(X)/n)))
    # G0 = ginv(X%*%t(X)/n)####different than matlab code: change inv to ginv
    Tn = rep(0, k_max)
    for (k in c(1:k_max)){
      sm = matrix(0, p, p)
      for (t in c(1:(n-k))){
        sm = sm + X[, t+k]%*%t(X[, t])
      }
      sm = sm/n
      Rmk = max(max(abs(sqrt(n)*R0 %*% sm %*% R0)))
      Tn[k] = Rmk
    }
    ###initial value of Tnn
    Tnn = rep(0, length(kk))
    for (j in c(1:length(kk))){
      Tnn[j] = max(Tn[1:kk[j]])
    }

    btemp = 1/sqrt(diag((X %*% t(X))/n))
    btemp = t(btemp)
    dm = c(1, rep(0, (p^2 -1)))

    for (j in c(1:p)){
      indp = c(((j-1)*p + 1):((j-1)*p + p))
      dm[indp] = btemp*btemp[j]
    }

    W = rep(dm, k_max)

    Sn = diag(rep(1, n), n, n)
    for (i in c(1:n)){
      for (j in c(1:n)){
        if (i != j){
          xu = (i-j)/bw
          Sn[i,j] = 25/(12*pi^2*xu^2)*(sin(6*pi*xu/5)/
                                         (6*pi*xu/5)-cos(6*pi*xu/5))
        }
      }
    }




    res = rep(0, length(kk))
    for (jj in c(1:length(kk))){
      K = kk[jj]
      nt = n - K

      ft = diag(1, p^2*K, nt)
      for (tt in c(1:nt)){
        for (j in c(1:K)){
          ind = ((j-1)*p^2+1):((j-1)*p^2+p^2)
          atemp = X[,tt+j] %*% t(X[,tt])
          ft[ind, tt] = rbind(atemp)
        }
      }
      mm = apply(ft, 1, mean)
      ft = ft - mm
      cv = aeS(ft = ft, Sn = Sn[1:nt, 1:nt], W = W[1: (p^2*K)], M, alpha = 0.05)
      res[jj] = (Tnn[jj]>cv)

    }


  }else if(type == 2){
    #test_LM
#     Y = par[[1]]
#     k = par[[2]]
    k = k_max
    p = dim(Y)[1]
    n = dim(Y)[2]
    Y = sweep(Y, MARGIN = 1, (apply(X = Y, MARGIN = 1, FUN = mean)), FUN = "-")
    YZz = c()
    for (j in (k:1)){
      YZz = rbind(YZz, Y[, j:(n-k+j)])
    }
    YZz = rbind(rep(1, n-k+1), YZz)
    nt = dim(YZz)[2]
    Yz = YZz[, 1:(nt-1)]
    Yy = Y[,(k+1):n]
    A =  (Yy %*% t(Yz)) %*% solve(Yz%*%t(Yz)) # A/B in matlab ==> B%*%solve(A)
    Z = Yy - A%*%Yz
    LM = (n - k)*(p - sum(diag(solve(Yy%*%t(Yy)) %*% (Z %*% t(Z))))) # A\B ==> solve(A)%*%B
    if (p > 25){
      Tstat1 = (LM - k*p^2)/sqrt(2*k*p^2)
      res = (abs(Tstat1) > qnorm(1- alpha/2))

    }else{
      res = (LM > qchisq(p = 1 - alpha, df = k*p^2))
    }
    # return (res)

  }else if(type == 3){
    #test_pre
    # Y = par[[1]]
    # k_max = par[[2]]
    # k = par[3]
    p = dim(Y)[1]
    n = dim(Y)[2]
    G0 = solve(Y%*%t(Y)/n)
    ll = rep(0, k_max)
    ls = rep(0, k_max)
    for (k in (1:k_max)){
      sm = matrix(0, p, p)
      for (t in (1:(n-k))){
        sm = sm + Y[, t+k] %*% t(Y[, t])
      }
      sm = sm/n
      ll[k] = sum(diag(t(sm) %*% t(G0) %*% sm %*% G0))
      ls[k] = ll[k]/(n-k)
    }
    q1 = q2 = q3 = rep(0, length(kk))
    q1n = q2n = q3n = rep(0, length(kk))
    for (ix in (1: length(kk))){
      j = kk[ix]
      Q1 = sum(ll[1:j])*n
      Q2 = sum(ls[1:j])*n^2
      Q3 = Q1 + p^2 * j * (j+1) / (2*n)
      cvK = qchisq(1-alpha, p^2*j)
      q1[ix] = Q1 > cvK
      q2[ix] = Q2 > cvK
      q3[ix] = Q3 > cvK

      Q1n = (Q1 - p^2*j)/sqrt(2*p^2*j)
      Q2n = (Q2 - p^2*j)/sqrt(2*p^2*j)
      Q3n = (Q3 - p^2*j)/sqrt(2*p^2*j)

      q1n[ix] = Q1n > qnorm(1-alpha);
      q2n[ix] = Q2n > qnorm(1-alpha);
      q3n[ix] = Q3n > qnorm(1-alpha);

    }
    res = cbind(q1,q2,q3, q1n, q2n, q3n)
  }else{
    #test_TB
    # Y = par[[1]]
    p = dim(Y)[1]
    n = dim(Y)[2]
    Y = sweep(Y, MARGIN = 1, apply(Y, 1, mean), FUN = "-")
    Y = rbind(t(rep(1, n)), Y)
    Yz = Y[,1:(n-1)]
    Yy = Y[, 2:n]
    A = (diag(rep(1, (p+1))) %*% Yy) %*% t(solve(Yz%*% t(Yz)) %*% Yz) # A\B ==> solve(A)%*%B

    Ytemp = A %*% Yz
    Z = Yy[2:(p+1), ] - Ytemp[2:(p+1),]

    U = max(det(X%*%t(X)/n) / det(Yy[2:(p+1),] %*% t(Yy[2:(p+1),])/n ), 0)
    Tstat = -log(U) * (n-p-3/2)
    if (p > 25){
      Tstat1 = (Tstat - p^2) / sqrt(2*p^2)
      res = (abs(Tstat1) > qnorm(1-alpha/2))

    }else{
      res = (Tstat > qchisq(1-alpha, p^2))
    }
  }

  return(res)
}

