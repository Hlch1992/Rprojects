modfit=function(datam,param){
  ########## get data
  yFdata=datam$yFdata;xFdata=datam$xFdata;conFdata=datam$conFdata;group=datam$group

  ########## estimated parameters
  beta=param$beta;gamma=param$gamma;sigma=param$sigma;cihatN2=param$cihatN2

  ########## estimated random parameters
  #;civarN2=ciest$civarN2

  ########## set dimentions
  N=length(yFdata);N2=length(unique(group));L=ncol(conFdata)


  mu <- exp(cbind(1,as.matrix(xFdata)) %*% beta+cihatN2[group])
  phi <- 1/(1+exp(-as.matrix(cbind(1,conFdata)) %*% gamma - cihatN2[(group+N2)]))

  ypred= (1-phi)*mu;yvar=(1-phi)*(phi*mu^2+mu)
  mse=mean((yFdata-ypred)^2)

  psse=unlist(sapply(1:N, function(x) ifelse(ypred[x]==0,0,(yFdata[x]-ypred[x])^2/yvar[x])))
  psmse=mean(psse)

  #if(sum(civarN2)!=0) {l=l-(cihatN2[group]^2/sigma[1]^2-log(2*pi*sigma[1]^2))/2-(cihatN2[(group+N2)]^2/sigma[2]^2-log(2*pi*sigma[2]^2))/2}
  #l[l==-Inf]=0

  return(list(rmse=mse,psmse=psmse,ypred=ypred))
}
