loglcal=function(mu,phi,yFdata){
  iszero=yFdata <= 0
  loglik0 <- log(phi + exp(log(1 - phi) - mu))
  loglik1 <- log(1 - phi) + dpois(yFdata, lambda = mu, log = TRUE)
  l=sum(sapply(1:length(iszero), function(l)  ifelse(iszero[l], loglik0[l],loglik1[l])))
}
