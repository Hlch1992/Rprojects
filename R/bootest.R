bootest=function(datalist=NA,paralist=NA,tunlist=NA){
  ########## get data
  yFdataNM=datalist$yFdataNM;xFdataNM=datalist$xFdataNM;conFdataNL=datalist$conFdataNL;group=datalist$group
  xFdataNM0=datalist$xFdataNM0

  ########## set dimentions
  N=nrow(yFdataNM);N2=length(unique(group));tmax=ceiling(N/N2);L=ncol(conFdataNL);M=ncol(yFdataNM)

  ########## estimated parameters
  beta=paralist$beta;gamma=paralist$gamma;sigma=paralist$sigma;ciest=paralist$ciest

  ########## fixed tunning parameter
  lambdabeta=tunlist$lambdabeta;lambdagamma=tunlist$lambdagamma;weight1all=tunlist$weight1all;weight2all=tunlist$weight2all
  selgamma=tunlist$selgamma

  betaboot=gammaboot=sigmaboot=ciestboot=NULL
  for(m in 1:M){
    datam=list(yFdata=yFdataNM[,m],xFdata=cbind(T0=xFdataNM0[,m],xFdataNM),conFdata=conFdataNL,group=group)
    param=list(betainit=beta[,m],gammainit=gamma[,m],sigmainit=sigma[m])
    weightm=list(w1=weight1all[,m],w2 = weight2all[,m])

    restemp = emnr(data = datam, para = param,family=family, weight = weightm,selgamma=selgamma)#,Lm=Lm

    betaboot=cbind(betaboot,restemp$para$beta);gammaboot=cbind(gammaboot,restemp$para$gamma);sigmaboot=c(sigmaboot,restemp$para$sigma)
    ciestboot=cbind(ciestboot,restemp$ciestm$cihatN2)
  }

  return(c(as.vector(betaboot),as.vector(gammaboot),sigmaboot))
}
