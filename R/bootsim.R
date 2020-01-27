bootsim=function(datalist=NA,paralist=NA){
  ########## get data
  yFdataNM=datalist$yFdataNM;xFdataNM=datalist$xFdataNM;conFdataNL=datalist$conFdataNL;group=datalist$group
  xFdataNM0=datalist$xFdataNM0

  ########## set dimentions
  N=nrow(yFdataNM);N2=length(unique(group));L=ncol(conFdataNL);M=ncol(yFdataNM)

  ########## estimated parameters
  beta=paralist$beta;gamma=paralist$gamma;sigma=paralist$sigma;biest=paralist$ciest

  logmu=as.matrix(cbind(1,xFdataNM)) %*% beta[-2,] + ciest[group] + t(t(xFdataNM0) * beta[2,])

  phi=1/(1+exp(-as.matrix(cbind(1,conFdataNL)) %*% gamma - ciest[(group+N2)]))

  phi[phi==.5]=0

  adlogmu=as.numeric(yFdataNM!=0) * logmu

  invtheta=colSums((yFdataNM-exp(adlogmu))^2-exp(adlogmu))/colSums(exp(adlogmu)^2)

  simyFdataNM=sapply(1:M, function(x) sapply(adlogmu[,x], function(t) rnbinom(1,(1/abs(invtheta[x])),mu=exp(t)))) * as.numeric(yFdataNM!=0)

  return(simyFdataNM)
}
