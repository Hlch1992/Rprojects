Update_p_M_old<-function(data,para,expectres,lambda1,lambda2,weight1,weight2,group,selgamma=F){

  ###### require packags
  suppressMessages(require(expm))

  ########## get data
  yFdata=data[[1]];xFdata=data[[2]];conFdata=data[[3]]

  ########## estimated parameters
  beta=para[[1]];gamma=para[[2]];sigma=para[[3]]

  ######### expectation in E-step
  bihatN2=expectres[[1]];bivarN2=expectres[[2]];uhatN=expectres[[3]]

  ########## set dimentions
  N=length(yFdata);N2=length(unique(group));tmax=ceiling(N/N2);L=ncol(conFdata)

  gammaold=gamma;betaold=beta;sigmaold=sigma

  if(sum(yFdata==0)==0) {gamma=rep(0,L+1)
  }else{
    if(ncol(as.matrix(conFdata))>1 & selgamma) {

      fitgamma=tryCatch({ withTimeout({
        glmnet(x=as.matrix(conFdata),y=cbind(1-uhatN,uhatN),family = 'binomial',lambda = lambda2,alpha = 1,penalty.factor = weight2)
      }, timeout = 60,onTimeout = 'error')}, error = function(e) {NA})

    }else{fitgamma=glm(y~.,family = binomial,data = data.frame(conFdata,y=uhatN))}
    gamma=as.vector(coef(fitgamma))
  }


  # Update_gamma(uhatN,conFdata,gammaold,lambda2,pen='lasso',weight = weight2)
  betasigma=Update_betasigma(data=data,beta=betaold,sigma=sigmaold,expectres=expectres,lambda1=lambda1,weight = weight1,group)
  beta=betasigma[[1]];sigma=betasigma[[2]]


  return(list(beta,gamma,sigma))
}
