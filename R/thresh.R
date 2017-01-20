thresh = function(Sigma_y, Y, my, k, n, p, delta){
  b = rep(0, p^2)
  s = 0
  p1 = p
  for (i in (1:p)){
    for (j in (1:p)){
      s = s+1
      theta = 0
      for (t in (1:(n-k))){
        theta = theta + ((Y[i, t+k] - my[i])*(Y[j,t] - my[j])
                         -Sigma_y[i,j])^2

      }
      theta = theta/n
      lam = delta*sqrt(theta*log(p1)/n)
      if (abs(Sigma_y[i,j]) < lam){
        Sigma_y[i,j] = 0
      }
      b[s] = Sigma_y[i,j]

    }
  }
  return(b)
}
