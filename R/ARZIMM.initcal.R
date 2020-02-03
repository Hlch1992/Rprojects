ARZIMM.initcal <-
  function(datalist,family='Poisson',iniw=TRUE,selectpara=list(selgamma=TRUE,selcri='BIC',pen='adalasso')){


  ##### read data
  yFdataNM=datalist$yFdataNM;xFdataNM=datalist$xFdataNM;conFdataNL=datalist$conFdataNL;group=datalist$group
  xFdataNM0=datalist$xFdataNM0

  ##### define basic parameters
  N2=length(unique(group));M=ncol(yFdataNM);N=nrow(xFdataNM);L=ncol(conFdataNL);q=ncol(xFdataNM)

  selgamma=selectpara$selgamma;selcri=selectpara$selcri;pen=selectpara$pen
  #################################################################
  ############### handling the initial value ######################
  #################################################################

  ######## define initial parameter matrix
  betainitMM=matrix(0,(q+2),M);gammainitLM = matrix(0,(L+1),M);sigmainitM = rep(0,2*M)

  ######## define weights
  weight1all=matrix(1,ncol = M,nrow = q+1) ;weight2all=matrix(1,ncol = M,nrow = L)
  for (m in 1:M) {

    yind=which(yFdataNM[,m]>0)

    ###### x data
    Xtemp=as.matrix(cbind(xFdataNM0[,m],xFdataNM)[yind,])

    ###### y data
    Ytemp=yFdataNM[yind,m]

    ###### lasso
    if(pen=='lasso') {
      fitini <- cv.glmnet(Xtemp,Ytemp,family = 'poisson')#nlambda=ntun,


      ###### adaptive lasso
    }else if(pen=='adalasso') {
      fit.ridge <- glmnet(Xtemp,Ytemp,family='poisson', alpha=0)#nlambda=ntun,
      # View(as.matrix(coef(fit.ridge)))
      BIC=log(length(Ytemp))*fit.ridge$df-2*apply(as.matrix(cbind(1,Xtemp)) %*% as.matrix(coef(fit.ridge)),2,function(ypred) sum(dpois(Ytemp,exp(ypred),log=TRUE)))
      w3 <- 1/abs(matrix(coef(fit.ridge,s=fit.ridge$lambda[which.min(BIC)])[-1, 1]))^1
      # w3
      w3[w3[,1] == Inf] <- w3[m+1] <- max(w3[w3!=Inf])*100
      w3[1]=0

      weight1all[,m]=w3

      fitini <- glmnet(Xtemp,Ytemp,family='poisson', alpha=1, penalty.factor=w3)
      # View(as.matrix(coef(fitini)))
    }

    ######## handleing tuning parameters
    BIC=log(length(Ytemp))*fitini$df-2*apply(as.matrix(cbind(1,Xtemp)) %*% as.matrix(coef(fitini)),2,function(ypred) sum(dpois(Ytemp,exp(ypred),log=TRUE)))
    BIC[which(as.matrix(coef(fitini))[1,] < 4 | as.matrix(coef(fitini))[1,] > 17)]=max(BIC)

    ######## passing the initial valie
    betainitMM[,m]=as.vector(coef(fitini,s=fitini$lambda[which.min(BIC)]))
    betainitMM[1,m]=mean(Xtemp[,(m+1)])*betainitMM[(m+2),m]+ betainitMM[1,m]
    betainitMM[(m+2),m]=0


    gmu_bi=log(yFdataNM[,m]+1)
    gmu_bi[yind]=cbind(1,Xtemp) %*% betainitMM[,m]

    sigmainitM[m]=sd(tapply((log(yFdataNM[,m]+1)-gmu_bi),group,mean))

    ###### concomitant variable initials
    if(sum(yFdataNM[,m]==0)>1){
      if(selgamma){
        ##### lasso
        if(pen=='lasso'){fit <- glmnet(as.matrix(conFdataNL), as.numeric(yFdataNM[,m]==0), family = "binomial")

        ##### adaptive lasso
        }else if(pen=='adalasso') {
          fit.ridge <- glmnet(as.matrix(conFdataNL), as.numeric(yFdataNM[,m]==0),family='binomial', alpha=0)
          BIC=log(N)*fit.ridge$df+2*apply(as.matrix(cbind(1,conFdataNL)) %*% as.matrix(coef(fit.ridge)),2,function(ypred)
            sum(log(1+exp(ypred))-as.numeric(yFdataNM[,m]==0)*ypred))
          w3 <- 1/abs(matrix(coef(fit.ridge, s=which.min(BIC))[-1, 1]))^1
          w3[w3[,1] == Inf] <- 999999999
          weight2all[,m]=w3
          fit <- glmnet(as.matrix(conFdataNL), as.numeric(yFdataNM[,m]==0),family='binomial',
                        alpha=1, penalty.factor=w3)
        }

        BIC=log(N)*fit$df+2*apply(as.matrix(cbind(1,conFdataNL)) %*% as.matrix(coef(fit)),2,function(ypred)
          sum(log(1+exp(ypred))-as.numeric(yFdataNM[,m]==0)*ypred))

        ######## initial values
        gammainitLM[,m]=as.vector(coef(fit,s=fit$lambda[which.min(BIC)]))


        gm1 <- tryCatch({
          glmer(as.formula(paste0('y ~ ',paste(colnames(conFdataNL)[which(gammainitLM[-1,m]!=0)],collapse = '+'),'+ (1 | group)')),
                data = data.frame(y=as.numeric(yFdataNM[,m]==0),conFdataNL,group=group),family = binomial,
                control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))

        }, error = function(e) {
          glmer(y~(1|group),data = data.frame(y=as.numeric(yFdataNM[,m]==0),conFdataNL,group=group),family = binomial,
                control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))

        })

        sigmainitM[(m+M)]=sqrt(as.numeric(VarCorr(gm1)))* gammainitLM[1,m]/coef(summary(gm1))[ 1, "Estimate"]

      }else{
        gm1 <- tryCatch({
          glmer(as.formula(paste0('y ~ ',paste(colnames(conFdataNL),collapse = '+'),'+ (1 | group)')),
                data = data.frame(y=as.numeric(yFdataNM[,m]==0),conFdataNL),family = binomial,
                control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
        }, error = function(e) {
          fit <- cv.glmnet(as.matrix(conFdataNL), as.numeric(yFdataNM[,m]==0),family='binomial', alpha=1)
          glmer(as.formula(paste0('y ~ ',paste(colnames(conFdataNL)[which(as.vector(coef(fit,s=fit$lambda.1se))[-1]!=0)],collapse = '+'),'+ (1 | group)')),
                data = data.frame(y=as.numeric(yFdataNM[,m]==0),conFdataNL,group=group),family = binomial,
                control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
        })


        gammainitLM[,m]=coef(summary(gm1))[ , "Estimate"]
        sigmainitM[(m+M)]=sqrt(as.numeric(VarCorr(gm1)))

      }

    }else if(sum(yFdataNM[,m]==0)==0){
      gammainitLM[,m]=0;sigmainitM[(m+M)]=0
    }else{
      gammainitLM[,m]=coef(glm(Y~.,family = binomial,data=data.frame(Y=as.numeric(yFdataNM[,m]==0),conFdataNL)))
      sigmainitM[(m+M)]=0
    }
  }

  gammainitLM[is.na(gammainitLM)]=0

  return(list(initpara=list(betaini=betainitMM,gammaini=gammainitLM,sigmaini=sigmainitM,weight1all=weight1all,weight2all=weight2all,iniw=iniw,calini=T),
              tunpara=list(lambdaseq1=NULL,lambdaseq2=NULL,ntun=50,epsilon = c(5e4,5))))
}

