sgTs_thre = function(X, k0, delta, opt=1, S1=1){
  n = dim(X)[1]
  p = dim(X)[2]
  if (opt == 1){
    M = cov(X)
    r = eigen(M)
    G = r$vectors
    ev = sqrt(r$values)
    D = matrix(rep(0, p^2), nrow = p)
    for (i in (1:p)){
      if (ev[i] > 0){
        D[i,i] = 1/(ev[i])
      }
    }
    M1 = (G) %*% D %*% t(G)
    X1 = M1 %*% t(X)
    Xn = X1
  }else if(opt== 2){
    M = S1
    X1 = solve(sqrtm(M)) %*% t(X) #A\B ==> solve(A)%*%B
    Xn = X1
    #   }else if(opt == 3){
    #     X = t(X)
    #     M =CLIME(X, sqrt(log(p)/n)*2)
    #     r = eigen(M)
    #     G = r$vectors
    #     ev = sqrt(r$values)
    #     M1 = G %*% diag(sqrt(ev), nrow = length(ev), ncol = length(ev)) %*% t(G)
    #     X1 = M1 %*% X
    #     Xn = X1
  }

  mean_Xn = apply(Xn, 1, mean)
  Wy = diag(1, p, p)

  if (p<6){
    for (k in (1:k0)){
      S = autocovm(t(Xn), k)
      Wy = Wy + (1 - k/(k0+1)) * S %*% t(S)
    }
  }else{
    for (k in (1:k0)){
      Sigma_y = autocovm(t(Xn), k)
      res = thresh(Sigma_y, Xn, mean_Xn, k, n, p, delta) #need to be done
      Sigma_ynew = matrix(res, p,p) %*% t(matrix(res,p,p))
      Wy = Wy + (1 - k/(k0 + 1)) * Sigma_ynew %*% t(Sigma_ynew)

    }
  }
  rr = eigen(Wy)
  G = rr$vectors
  ev = rr$values
  X1 = t(t(G) %*%Xn)
  return(X1)

}
