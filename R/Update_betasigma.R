Update_betasigma<-function(data,beta,sigma,expectres,lambda1,weight,group){
  ###### require packags
  suppressMessages(require(expm))

  ########## get data
  yFdata=data[[1]];xFdata=data[[2]];conFdata=data[[3]]

  ######### expectation in E-step
  bihatN2=expectres[[1]];bivarN2=expectres[[2]];uhatN=expectres[[3]]

  ########## set dimentions
  N=length(yFdata);N2=length(unique(group));tmax=ceiling(N/N2);L=ncol(conFdata)

  sigma=sd(bihatN2)
  ##########
  eps=10
  iter=1

  while(eps>10^(-3) & iter<=10){

    betaold=beta#;sigmaold=sigma
    linkmuN=cbind(1,as.matrix(xFdata)) %*% beta + as.vector(bihatN2[group])
    ###### negative second derivative of log likelihood: -l''(bi)
    diagsDbiN=(1-uhatN)*exp(linkmuN)
    diagsDbiN[is.na(diagsDbiN)]=10^7
    diagsDbiN[diagsDbiN==Inf]=10^7
    risigma=tapply(diagsDbiN, group, sum)
    ri=risigma+1/sigma^2

    fDbeta=as.vector(diagsDbiN*(1+1/as.vector(ri[group])/2)-(1-uhatN)*yFdata) %*% cbind(1,as.matrix(xFdata))


    #firstDbeta(xFdata,yFdata,betaold,sigmaold,uhatN,ri,linkmuN)

    wtsdbeta=as.vector(diagsDbiN*(1+1/2/as.vector(ri[group])-as.vector((risigma/ri^2)[group])/2))

    # wtsdbeta=as.vector((1-uhatN)*exp(linkmuN)*(1+1/2/rep(ri,tmax)-rep(risigma/ri^2,tmax)/2))



    sDbeta=matrix(0,ncol(xFdata)+1,ncol(xFdata)+1)
    for (i in 1:N) {
      sDbeta=sDbeta+wtsdbeta[i]*outer(c(1,unlist(xFdata[i,])),c(1,unlist(xFdata[i,])))
    }


    # secondDbeta(xFdata,yFdata,betaold,sigmaold,uhatN,negsdl,linkmuN)
    if(sum(sDbeta)==0) {break}
    sDroot=tryCatch({sqrtm(sDbeta)}, error = function(e) {sqrtsingular(sDbeta)})

    if(is.complex(sDroot)) {if(max(abs(Im(sDroot)))<1e4) sDroot=Re(sDroot)}

    my_Y=tryCatch({solve(sDroot+10^(-6)*diag(length(beta)))%*%t((as.vector(sDbeta%*%betaold)-fDbeta)/sqrt(2))
    },error = function(e){solve(sDroot+min(sDroot)*diag(length(beta)))%*%t((as.vector(sDbeta%*%betaold)-fDbeta)/sqrt(2))
    })
    my_X=sDroot/sqrt(2)
    fit=glmnet(my_X,my_Y,family="gaussian",alpha=1,lambda=lambda1,penalty.factor = weight,intercept=FALSE)
    beta=c(as.matrix(fit$beta))


    eps=sum(abs(c(betaold-beta)))#,sigmaold-sigma
    #  print(paste0('updatbetasigma',eps))
    iter=iter+1
    if(max(abs(beta))>20 |  max(ri)>1e9 | eps>20) {beta=betaold;break}#sigma>20 |
    #;sigmaold=sigma

  }

  return(list(beta,sigma))
}

