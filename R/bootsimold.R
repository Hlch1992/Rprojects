bootsimold=function(datalist=NA,paralist=NA){
  ########## get data
  yFdataNM=datalist$yFdataNM;xFdataNM=datalist$xFdataNM;conFdataNL=datalist$conFdataNL;group=datalist$group

  ########## set dimentions
  N=nrow(yFdataNM);N2=length(unique(group));tmax=ceiling(N/N2);L=ncol(conFdataNL)

  ########## estimated parameters
  beta=paralist$beta;gamma=paralist$gamma;sigma=paralist$sigma;biest=paralist$biest

  logmu=as.matrix(cbind(1,xFdataNM)) %*% beta + biest[group,]

  phi=1/(1+exp(-as.matrix(cbind(1,conFdataNL)) %*% gamma))

  phi[phi==.5]=0

  adlogmu=as.numeric(yFdataNM!=0) * logmu

  simyFdataNM=apply(adlogmu, 2, function(x) unlist(sapply(x, function(t) rpois(1,exp(t))))) * as.numeric(yFdataNM!=0)

  return(simyFdataNM)
}
