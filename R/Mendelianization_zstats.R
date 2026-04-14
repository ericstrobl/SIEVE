Mendelianization_zstats <-function(Alpha,Sigma,Zstat,leads){
  
  d = ncol(as.matrix(Alpha))
  Z = c()
  for (i in 1:length(leads)) {
    if (d > 1){
      a <- Alpha[, leads[i], drop = FALSE]      # alpha_j (m x 1)
    } else{
      a <- Alpha
    }
    denom_i <- sqrt(pmax(as.numeric(t(a) %*% Sigma %*% a), 1e-20))
    
    Z <- cbind(Z,as.vector((Zstat %*% a) / denom_i))
    
  }
  p = 2 * pnorm(abs(Z), lower.tail = FALSE)
  
  return(list(Z=Z,p=p))
}
