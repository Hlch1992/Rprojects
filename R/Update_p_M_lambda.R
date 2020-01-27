Update_p_M_lambda<-function(data,para,expectres,lambda,weight,selgamma=F,sel=F){

  ###### require packags
  suppressMessages(require(expm))

  ########## get data
  yFdata=data$yFdata;xFdata=data$xFdata;conFdata=data$conFdata;group=data$group

  ########## estimated parameters
  beta=para$beta;gamma=para$gamma;sigma=para$sigma

  ######### expectation in E-step
  cihatN2=expectres$cihatN2;civarN2=expectres$civar;uhatN=expectres$uhatN

  ########## set dimentions
  N=length(yFdata);N2=length(unique(group));L=ncol(conFdata)

  gammaold=gamma;betaold=beta;sigmaold=sigma

  if(sum(yFdata==0)==0) {
    gamma=rep(0,L+1);sigma[2]=0
    sigma[1]=sd(cihatN2[1:N2])

    fitbeta=tryCatch({ withTimeout({
      glmnet(as.matrix(xFdata),yFdata,family = 'poisson',penalty.factor = weight$w1,lambda = lambda[1],
             offset = cihatN2[group],alpha = 1)
    }, timeout = 60,onTimeout = 'error')}, error = function(e) {NA})

    beta=as.vector(coef(fitbeta))

  }else{
    sigma[2]=sd(cihatN2[(which(cihatN2[(1+N2):(2*N2)]!=0)+N2)])
    if(ncol(as.matrix(conFdata))>1 & selgamma) {
      fitgamma=tryCatch({ withTimeout({
        glmnet(x=as.matrix(conFdata),y=cbind(1-uhatN,uhatN),family = 'binomial',alpha = 1,lambda = lambda[2],
               offset = cihatN2[(group+N2)],penalty.factor = weight$w2)#
      }, timeout = 60,onTimeout = 'error')}, error = function(e) {NA})

      gamma=as.vector(coef(fitgamma))
    }else{fitgamma=glm(y~.+offset(cihatN2[(group+N2)]),family = binomial,data = data.frame(conFdata,y=uhatN))
    gamma=coef(fitgamma)}


    if(sum(cihatN2[1:N2]!=0)<=1) {
      sigma[1]=0

      fitbeta=tryCatch({ withTimeout({
        glmnet(as.matrix(xFdata[which(uhatN!=1),]),yFdata[which(uhatN!=1)],family = 'poisson',lambda = lambda[1],
               penalty.factor = weight$w1,alpha = 1)
      }, timeout = 60,onTimeout = 'error')}, error = function(e) {NA})


      beta=as.vector(coef(fitbeta))

    }else{
      sigma[1]=sd(cihatN2[which(cihatN2[1:N2]!=0)])
      fitbeta=tryCatch({ withTimeout({
        glmnet(as.matrix(xFdata[which(uhatN!=1),]),yFdata[which(uhatN!=1)],family = 'poisson',lambda = lambda[1],
               penalty.factor = weight$w1,offset = (cihatN2[group])[which(uhatN!=1)],alpha = 1)
      }, timeout = 60,onTimeout = 'error')}, error = function(e) {NA})




      beta=as.vector(coef(fitbeta))
    }
  }
  if(sum(is.na(sigma))>0) sigma[is.na(sigma)]=0

  return(list(beta=beta,gamma=gamma,sigma=sigma))
}
