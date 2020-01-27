loglall=function(datalist,paralist){
  ########## get data
  yFdata=datalist$yFdataNM;xFdata=datalist$xFdataNM;conFdata=datalist$conFdataNL;group=datalist$group

  ########## estimated parameters
  beta=paralist$beta;gamma=paralist$gamma;sigma=paralist$sigma;cihatN2=paralist$ciest

  ########## estimated random parameters
  #;civarN2=ciest$civarN2
  ########## set dimentions
  N=nrow(yFdata);N2=length(unique(group));L=ncol(conFdata);M=ncol(yFdata)
  ypred=l=matrix(0,N,M);mse=psmse=rep(0,M)

  mu <- exp(cbind(1,as.matrix(xFdata)) %*% beta+cihatN2[group,])
  phi <- 1/(1+exp(-as.matrix(cbind(1,conFdata)) %*% gamma - cihatN2[(group+N2),]))

  iszero=yFdata <= 0
  loglik0 <- log(phi + exp(log(1 - phi) - mu))
  loglik1 <- log(1 - phi) + sapply(1:M, function(indcol) dpois(yFdata[,indcol], lambda = mu[,indcol], log = TRUE))

  l=outer(1:dim(iszero)[1],1:dim(iszero)[2], function(x,y) {
    foo=function(x,y) ifelse(iszero[x,y], loglik0[x,y],loglik1[x,y])
    mapply(foo, x,y)
  })
  ypred= (1-phi)*mu
  return(list(nl=l,ypred=ypred))
}
