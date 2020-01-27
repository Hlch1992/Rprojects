glmaicbic=function(fit) {
  tLL <- fit$nulldev - deviance(fit)
  k <- fit$df
  n <- fit$nobs
  AIC <- -tLL+2*k+2*k
  BIC<-log(n)*k - tLL
  return(rbind(BIC=BIC,AIC=AIC))
}
