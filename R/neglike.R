neglike=function(datam,param){
  ########## get data
  yFdata=datam$yFdata;xFdata=datam$xFdata;conFdata=datam$conFdata;group=datam$group

  ########## estimated parameters
  beta=param$beta;gamma=param$gamma;sigma=param$sigma;cihatN2=param$cihatN2

  ########## estimated random parameters
  #;civarN2=ciest$civarN2

  ########## set dimentions
  N=length(yFdata);N2=length(unique(group));L=ncol(conFdata)


  mu <- exp(as.vector(cbind(1,as.matrix(xFdata)) %*% beta)+cihatN2[group])
  phi <- as.vector(1/(1+exp(-as.matrix(cbind(1,conFdata)) %*% gamma - cihatN2[(group+N2)])))

  iszero=yFdata <= 0
  loglik0 <- log(phi + exp(log(1 - phi) - mu))
  loglik1 <- log(1 - phi) + dpois(yFdata, lambda = mu, log = TRUE)
  l=sapply(1:length(iszero), function(l)  ifelse(iszero[l], loglik0[l],loglik1[l]))

  ypred= (1-phi)*mu
  mse=sqrt(mean((yFdata-ypred)^2))

  psse=unlist(sapply(1:N, function(x) ifelse(ypred[x]==0,0,(yFdata[x]-ypred[x])^2/ypred[x])))
  psmse=sqrt(mean(psse))

  #if(sum(civarN2)!=0) {l=l-(cihatN2[group]^2/sigma[1]^2-log(2*pi*sigma[1]^2))/2-(cihatN2[(group+N2)]^2/sigma[2]^2-log(2*pi*sigma[2]^2))/2}
  #l[l==-Inf]=0

  return(list(nl=l,rmse=mse,psmse=psmse,ypred=ypred))
}
