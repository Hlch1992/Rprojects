Update_p_M<-function(data,para,expectres,weight,selgamma=FALSE,sel=FALSE,fast=TRUE){

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
    if(!fast){
      fit.ridge <- tryCatch({ withTimeout({
        glmnet(as.matrix(xFdata),yFdata,family='poisson', offset = cihatN2[group],alpha=0)#nlambda=ntun,
      }, timeout = 60,onTimeout = 'error')}, error = function(e) {NA})
      # View(as.matrix(coef(fit.ridge)))
      BIC=log(length(yFdata))*fit.ridge$df-2*apply(as.matrix(cbind(1,xFdata)) %*% as.matrix(coef(fit.ridge)) + cihatN2[group],2,function(ypred)
        sum(dpois(yFdata,exp(ypred),log=TRUE)))
      weight$w1 <- 1/abs(matrix(coef(fit.ridge,s=fit.ridge$lambda[which.min(BIC)])[-1, 1]))^1
      # w3
      weight$w1[weight$w1[,1] == Inf] <- max(weight$w1[weight$w1!=Inf])*100

    }

    fitbeta=tryCatch({ withTimeout({
      glmnet(as.matrix(xFdata),yFdata,family = 'poisson',penalty.factor = weight$w1,#lambda = lambda[1],
             offset = cihatN2[group],alpha = 1)
    }, timeout = 60,onTimeout = 'error')}, error = function(e) {NA})


    if(is.na(fitbeta)) {return(list(beta=NA,gamma=NA,sigma=NA,lambda=NA,logl=NA))}

    betatmp=tryCatch({
      as.matrix(coef(fitbeta))
    }, error = function(e) {
      NA
    })

    if(is.na(betatmp)) {return(list(beta=NA,gamma=NA,sigma=NA,lambda=NA,logl=NA))}
    mu <- exp(cbind(1,as.matrix(xFdata)) %*% betatmp +cihatN2[group])
    phi <- 1/(1+exp(-as.matrix(cbind(1,conFdata)) %*% gamma - cihatN2[(group+N2)]))

    logltun=sapply(1:ncol(mu),function(x) loglcal(mu[,x],phi,yFdata))
    if(sel) {
      BIC=sapply(1:ncol(mu),function(x) log(N)*(fitbeta$df[x])-2*logltun[x])
      indtmp=which.min(BIC)
    }else{

      indtmp=which.max(logltun)
    }


  }else{
    sigma[2]=sd(cihatN2[(which(cihatN2[(1+N2):(2*N2)]!=0)+N2)])
    if(ncol(as.matrix(conFdata))>1 & selgamma) {
      if(!fast) {
        fit.ridge <- tryCatch({ withTimeout({
          glmnet(x=as.matrix(conFdata),y=cbind(1-uhatN,uhatN),family='binomial', offset = cihatN2[(group+N2)], alpha=0)
        }, timeout = 60,onTimeout = 'error')}, error = function(e) {NA})


        BIC=log(N)*fit.ridge$df+2*apply(as.matrix(cbind(1,conFdata)) %*% as.matrix(coef(fit.ridge)) + cihatN2[(group+N2)],2,function(ypred)
          sum(log(1+exp(ypred))-uhatN*ypred))
        weight$w2 <- 1/abs(matrix(coef(fit.ridge, s=which.min(BIC))[-1, 1]))^1
        weight$w2[weight$w2[,1] == Inf] <- 999999999
      }
      fitgamma=tryCatch({ withTimeout({
        cv.glmnet(x=as.matrix(conFdata),y=cbind(1-uhatN,uhatN),family = 'binomial',alpha = 1,#lambda = lambda[2],
                  offset = cihatN2[(group+N2)],penalty.factor = weight$w2)
      }, timeout = 60,onTimeout = 'error')}, error = function(e) {NA})


      if(is.na(fitgamma)) {return(list(beta=NA,gamma=NA,sigma=NA,lambda=NA,logl=NA))}
      gamma=as.matrix(coef(fitgamma,s=fitgamma$lambda.1se))

    }else{fitgamma=glm(y~.+offset(cihatN2[(group+N2)]),family = binomial,data = data.frame(conFdata,y=uhatN))
    gamma=coef(fitgamma)}


    if(sum(cihatN2[1:N2]!=0)<=1) {
      sigma[1]=0
      if(!fast){
        fit.ridge <- tryCatch({ withTimeout({
          glmnet(as.matrix(xFdata[which(uhatN!=1),]),yFdata[which(uhatN!=1)],family='poisson',alpha=0)#nlambda=ntun,
        }, timeout = 60,onTimeout = 'error')}, error = function(e) {NA})

        # View(as.matrix(coef(fit.ridge)))
        BIC=log(sum(uhatN!=1))*fit.ridge$df-2*apply(as.matrix(cbind(1,xFdata[which(uhatN!=1),])) %*% as.matrix(coef(fit.ridge)) ,2,function(ypred)
          sum(dpois(yFdata[which(uhatN!=1)],exp(ypred),log=TRUE)))
        weight$w1 <- 1/abs(matrix(coef(fit.ridge,s=fit.ridge$lambda[which.min(BIC)])[-1, 1]))^1
        # w3
        weight$w1[weight$w1[,1] == Inf] <- max(weight$w1[weight$w1!=Inf])*100

      }

      fitbeta=tryCatch({ withTimeout({
        glmnet(as.matrix(xFdata[which(uhatN!=1),]),yFdata[which(uhatN!=1)],family = 'poisson',
               penalty.factor = weight$w1,alpha = 1)
      }, timeout = 60,onTimeout = 'error')}, error = function(e) {NA})


      if(is.na(fitbeta)) {return(list(beta=NA,gamma=NA,sigma=NA,lambda=NA,logl=NA))}

      betatmp=tryCatch({
        as.matrix(coef(fitbeta))
      }, error = function(e) {
        NA
      })

      if(is.na(betatmp)) {return(list(beta=NA,gamma=NA,sigma=NA,lambda=NA,logl=NA))}

    }else{
      sigma[1]=sd(cihatN2[which(cihatN2[1:N2]!=0)])
      if(!fast){
        fit.ridge <- tryCatch({ withTimeout({
          glmnet(as.matrix(xFdata[which(uhatN!=1),]),yFdata[which(uhatN!=1)],family='poisson',
                 offset = (cihatN2[group])[which(uhatN!=1)],alpha=0)#nlambda=ntun,
        }, timeout = 60,onTimeout = 'error')}, error = function(e) {NA})



        # View(as.matrix(coef(fit.ridge)))
        BIC=log(sum(uhatN!=1))*fit.ridge$df-2*apply(as.matrix(cbind(1,xFdata[which(uhatN!=1),])) %*% as.matrix(coef(fit.ridge)) +
                                                      (cihatN2[group])[which(uhatN!=1)] ,2,function(ypred)
                                                        sum(dpois(yFdata[which(uhatN!=1)],exp(ypred),log=TRUE)))
        weight$w1 <- 1/abs(matrix(coef(fit.ridge,s=fit.ridge$lambda[which.min(BIC)])[-1, 1]))^1
        # w3
        weight$w1[weight$w1[,1] == Inf] <- max(weight$w1[weight$w1!=Inf])*100

      }

      fitbeta=tryCatch({ withTimeout({
        glmnet(as.matrix(xFdata[which(uhatN!=1),]),yFdata[which(uhatN!=1)],family = 'poisson',#lambda = lambda[1],
               penalty.factor = weight$w1,offset = (cihatN2[group])[which(uhatN!=1)],alpha = 1)
      }, timeout = 60,onTimeout = 'error')}, error = function(e) {NA})


      if(is.na(fitbeta)) {return(list(beta=NA,gamma=NA,sigma=NA,lambda=NA,logl=NA))}

      betatmp=tryCatch({
        as.matrix(coef(fitbeta))
      }, error = function(e) {
        NA
      })

      if(is.na(betatmp)) {return(list(beta=NA,gamma=NA,sigma=NA,lambda=NA,logl=NA))}

    }



    mu <- exp(cbind(1,as.matrix(xFdata)) %*% betatmp+cihatN2[group])
    phi <- 1/(1+exp(-as.matrix(cbind(1,conFdata)) %*% gamma - cihatN2[(group+N2)]))

    logltun=sapply(1:ncol(mu),function(x) loglcal(mu[,x],phi,yFdata))
    if(sel) {
      BIC=sapply(1:ncol(mu),function(x) log(N)*(fitbeta$df[x])-2*logltun[x])
      indtmp=which.min(BIC)
    }else{

      indtmp=which.max(logltun)
    }
  }




  if(sum(is.na(sigma))>0) sigma[is.na(sigma)]=0

  if(sum(yFdata==0)==0) {
    return(list(beta=betatmp[,indtmp],gamma=gamma,sigma=sigma,lambda=c(fitbeta$lambda[indtmp],0),logl=logltun[indtmp]))
  }else{
    return(list(beta=betatmp[,indtmp[1]],gamma=gamma,sigma=sigma,
                lambda=c(fitbeta$lambda[indtmp],fitgamma$lambda.1se),logl=logltun[indtmp]))}
}
