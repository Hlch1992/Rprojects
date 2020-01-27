emnrlambda=function(data,para,lambda,weight,family,selgamma=F){

  ########## get data
  yFdata=data$yFdata;xFdata=data$xFdata;conFdata=data$conFdata;group=data$group

  ########## set dimentions
  N=length(yFdata);N2=length(unique(group));L=ncol(conFdata)

  ########## estimated parameters
  beta=para$betainit;gamma=para$gammainit;sigma=para$sigmainit
  if(sum(yFdata==0)==0) gamma=rep(0,L+1)

  ########## set boundary
  BIG=1.0e2; SMALL=1.0e-2; converge=0;maxits=1000;its2=0; fail=0;

  ##########
  cihatN2=rep(0,2*N2)

  lnew=-Inf

  ########## start iterate
  while(its2 < maxits) {
    its2=its2+1
    betaold=beta; gammaold=gamma; sigmaold=sigma;cihatN2old=cihatN2;lold=lnew

    ###### score before transform for possion and logistic part
    theta_bi=cbind(1,as.matrix(xFdata)) %*% beta; logitp_ai=as.matrix(cbind(1,conFdata)) %*% gamma

    ###### E-step:
    ###### expectation of U ; bi ; var(bi)
    expectres=Updata_E(yFdata,theta_bi,logitp_ai,sigma,cihatN2,group)
    if (expectres$conv==1) {cihatN2=expectres$cihatN2;civarN2=expectres$civar;uhatN=expectres$uhatN}

    lnew=loglcal(mu=theta_bi+cihatN2[group],phi = 1/(1+exp(-logitp_ai-cihatN2[(group+N2)])),yFdata = yFdata)
    # print(lnew)
    ###### M-step:
    maximres=Update_p_M(data=data,para=list(beta=beta,gamma=gamma,sigma=sigma),#lambda=lambda,
                        expectres=expectres,weight=weight,selgamma=selgamma)
    beta=maximres$beta;gamma=maximres$gamma;sigma=maximres$sigma
    if(sum(is.na(c(beta,gamma,sigma)))>0) {converge=2;break}

    ####### test convergency
    eps=sum(abs(c(beta- betaold,gamma- gammaold,sigma- sigmaold)))#
    if(its2==1) epsold=eps
    if (eps<SMALL) {converge=1;break
    }else if(eps>5*BIG | eps>epsold*10 | max(abs(sigma))>BIG | min(sigma) < 0 | abs(beta[1]) > 50 | sum(beta)==0 | max(abs(cihatN2))> 50 | lnew < lold) {
      converge=2;break}
    epsold=eps
  }

  if(converge==2) {   beta = betaold;gamma = gammaold;sigma = sigmaold; cihatN2=cihatN2old}
  ###### score before transform for possion and logistic part
  theta_bi=cbind(1,as.matrix(xFdata)) %*% beta; logitp_ai=as.matrix(cbind(1,conFdata)) %*% gamma

  ###### E-step:
  ###### expectation of U ; bi ; var(bi)
  expectres=Updata_E(yFdata,theta_bi,logitp_ai,sigma,cihatN2,group)
  cihatN2=expectres$cihatN2;civarN2=expectres$civar;uhatN=expectres$uhatN
  maximres=Update_p_M_lambda(data=data,para=list(beta=beta,gamma=gamma,sigma=sigma),lambda=lambda,
                             expectres=expectres,weight,selgamma=selgamma)
  beta=maximres$beta;gamma=maximres$gamma;sigma=maximres$sigma


  para=list(beta=beta,gamma=gamma,sigma=sigma,cihatN2=cihatN2)
  mf=modfit(data,para)
  rmse=mf$rmse;psmse=mf$psmse;ypred=mf$ypred
  return(list(para=para,conv=converge,df=sum(unlist(para)!=0),mse=c(rmse,psmse),ypred=ypred))
}
