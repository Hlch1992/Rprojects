myparaset=function(tmp,wpw){
  library(stringr)
  
  
  nsimf=tmp[7];nsimt=tmp[8]
  
  ##### number of individuals
  N2=tmp[1] # 20; 50
  ##### number of time points
  timeN=tmp[2] # 10; 20
  ##### number of Microbiome
  M=tmp[3] # 20; 50
  ##### has mix effect or not
  mixeff=tmp[4] #1; 0
  ##### zero inflation proportion
  zeropro=tmp[5] #0.3; 0.6; 0
  effss=tmp[6]
  ##### informative interaction proportion
  ineffpro=.05
  
  #if(tmp[6]==0.01) {effss=0.1;effs=0.1}else if(tmp[6]==0.1) {effss=0.1;effs=1}else if(tmp[6]==1) {effss=1;effs=1}
  
  ##### set seed
  if(M==20) {seed=2010200}
  
  
  if(zeropro==0) {
    set.seed(seed)
    ########## baseline
    ########## Microbiome log mean
    interceptM=rnorm(M,9,2)
    
    ########## Microbiome random effect sigma
    if(mixeff==1) {sigmaM=10^((interceptM-9)/6-2)} else {sigmaM=rep(0,M)}
    ########## individual random intercept
    biN2M=sapply(sigmaM, function(x) rnorm(N2,0,x))
    
    ########## Individual microbiome counts
    baseFdataN2M=outer(1:N2,1:M,function(x,y) {
      foo=function(x,y) rpois(1,exp(interceptM[y]+biN2M[x,y]))
      mapply(foo, x, y)
    })
    rownames(baseFdataN2M)=paste0('subj',str_pad(1:N2, 2, pad = "0"))
    colnames(baseFdataN2M)=paste0('M',c(1:M))
    
    ########## zero inflation 
    ########## concomitant variable
    ind=rbinom(M*3,1,.5)
    conFdataN2L=cbind(W1=rnorm(N2,0,1),W2=rnorm(N2,0,1),W3=rbinom(N2,1,.33),W4=rbinom(N2,1,.33),
                      W5=floor(runif(N2,1,4)),W6=floor(runif(N2,1,4)))
    gammaLM=NA
    
    ########## interaction matrix
    betaMM=rep(0,M*M)
    
    if(mixeff==0){ 
      intind=c(32 , 93 ,369,  180 ,362 ,
               92, 259, 89 ,165 , 33, 
               248 , 99, 163 , 9 ,241, 
               233 ,221 ,363 , 17, 172 )
      effectsize=sqrt(effss)*c(0.8621417, 0.608599001, 0.22568930, 0.46165429, 0.31791238, 
                               0.22191497, 0.55076096, 0.50482682, 0.36285228, 1.0530587,
                               0.26325326, 0.22950178, 0.43863914, 1.40365321, 0.50312829, 
                               0.15789547, 0.23545650, 0.17002869, 1.20068537, 0.34281308)
      direction=c( 1,  1,  1, 1,  -1, -1 ,-1 ,-1 ,-1 ,-1 ,1 , 1 , 1, -1,  1,  1  ,-1  ,1, 1,  -1)
      
      
      
    }else if(mixeff==1 ){ 
      intind=c(32 , 100 ,378,  172 ,371 ,
               92, 258, 89 ,168 , 33, 
               249 , 97, 163 , 8 ,241, 
               238 ,221 ,374 , 17, 179 )
      effectsize=sqrt(effss)*c(0.8621417, 0.588599001, 0.22568930, 0.50165429, 0.31791238, 
                               0.26191497, 0.55076096, 0.50482682, 0.34285228, 1.0530587,
                               0.26325326, 0.22950178, 0.43863914, 1.40365321, 0.50312829, 
                               0.15789547, 0.23545650, 0.17002869, 1.20068537, 0.34281308)
      direction=c( 1,  1,  1, 1,  -1, -1 ,-1 ,-1 ,-1 ,-1 ,1 , 1 , 1, -1,  1,  1  ,-1  ,1, 1,  -1)
      
    }
    
    betaMM[intind]=effectsize*direction
    betaMM=matrix(betaMM,ncol = M,dimnames = list(paste0('M',c(1:M)),paste0('M',c(1:M))))

    
  }else if (zeropro==1){
    #setwd(wpw)
    if(N2==20) {load(paste0(wpw,"/Data/BLsub20m20.RData"))
    }else if(N2==50) {load(paste0(wpw,"/Data/BLsub50m20.RData"))}
    
    interceptM=interceptM-2
    set.seed(seed) 
    ########## Microbiome random effect sigma
    if(mixeff==1) {sigmaM=10^((interceptM-9)/6-2)} else {sigmaM=rep(0,M)}
    
    ########## individual random intercept
    biN2M=sapply(sigmaM, function(x) rnorm(N2,0,x))
    
    ########## zero inflation
    ziprop=1/(1+exp(-(cbind(1,conFdataN2L) %*% gammaLM)))
    ziprop[,which(colSums(gammaLM)==0)]=0
    uN2M=outer(1:N2,1:M,function(x,y) {
      foo=function(x,y) rbinom(1,1,ziprop[x,y])
      mapply(foo, x,y)
    })
    
    baseFdataN2M=outer(1:N2,1:M,function(x,y) {
      foo=function(x,y) rpois(1,exp(interceptM[y]+biN2M[x,y]))*(1-uN2M[x,y])
      mapply(foo, x, y)
    })
    rownames(baseFdataN2M)=paste0('subj',str_pad(1:N2, 2, pad = "0"))
    colnames(baseFdataN2M)=paste0('M',c(1:M))
    
    
    ########## interaction matrix
    betaMM=rep(0,M*M)
    
    if(mixeff==0){# & effss==1
      intind=c(222, 339,  40, 140, 341, 
               250, 102,  38, 314, 120,
               111, 138, 310, 230, 330,
               396, 365,  20,  10, 382)#sample(1:M^2,ineffpro*M*M,replace = F)
      effectsize=sqrt(effss)*c(0.188235605, 0.234132241, 0.41366703, 0.61070656, 0.06591708,
                                    0.04580601, 0.82368206, 0.35504387, 1.23225201, 0.97013089,
                                    0.06018784, 0.52154596, 0.99354102, 0.14859757, 0.105828898,
                                    0.37905519, 0.08377116, 0.88736514, 0.70141201, 0.31951335)#10^rnorm(ineffpro*M*M,-.5,.5)
      direction=c(1,  1,  1,  1,  1,
                  1, -1, -1,  1,  1,
                  -1, -1, -1, -1, -1,
                  -1, -1, -1,  1,  1)
      
    
    }else if(mixeff==1) {    # & effss==1
      
      intind=c(223, 339,  39, 140, 342, 
               250, 102,  37, 314, 120,
               111, 138, 320, 230, 332,
               396, 370,  20,  10, 382)#sample(1:M^2,ineffpro*M*M,replace = F)
      effectsize=sqrt(effss)*c(0.188235605, 0.234132241, 0.41366703, 0.61070656, 0.10591708,
                               0.04580601, 0.82368206, 0.35504387, 1.23225201, 0.90013089,
                               0.06018784, 0.52154596, 0.99354102, 0.14859757, 0.105828898,
                               0.37905519, 0.08377116, 0.88736514, 0.70141201, 0.31951335)#10^rnorm(ineffpro*M*M,-.5,.5)
      direction=c(1,  1,  1,  1,  1,
                  1, -1, -1,  1,  1,
                  -1, -1, -1, -1, -1,
                  -1, -1, -1,  1,  1)#sample(rep(c(1,-1),ineffpro*M*M/2),replace = F)
      
      }
      
    
    betaMM[intind]=effectsize*direction
    betaMM=matrix(betaMM,ncol = M,dimnames = list(paste0('M',c(1:M)),paste0('M',c(1:M))))
    
  }
  
  
  
  return(paralist=list(baseFdataN2M=baseFdataN2M,conFdataN2L=conFdataN2L,timeN=timeN,
                       interceptM=interceptM,betaMM=betaMM,gammaLM=gammaLM,sigmaM=sigmaM,biN2M=biN2M))
}



initcal=function(tmp,Varname=NA,Conname=NA,fdata=NA,IDname='ID',Tname='Time',family='Poisson',pen='adalasso',iniw=T,group=NA,
                 tunpara=list(lambdaseq='selfdefine',lambdaseq1=NULL,lambdaseq2=NULL,ntun=5,epsilon = c(1e6,1e4)),
                 selectpara=list(selgamma=T,selcri='BIC',adjcov=FALSE,Covname=NA)){

  ##### number of individuals
  N2=tmp[1] # 20; 50
  ##### number of time points
  timeN=tmp[2] # 10; 20
  ##### number of Microbiome
  M=tmp[3] # 20; 50
  ##### has mix effect or not
  mixeff=tmp[4] #1; 0
  ##### zero inflation proportion
  zeropro=tmp[5] #1; 0
  ##### effect size scale
  effss=tmp[6]
  ##### informative interaction proportion
  ineffpro=.05
  
  
  if(TRUE) {
    #betaini=initpara$betaini;gammaini=initpara$gammaini;sigmaini=initpara$sigmaini
    #weight1all=initpara$weight1all;weight2all=initpara$weight2all;iniw=initpara$iniw
    
    lambdaseq=tunpara$lambdaseq;lambdaseq1=tunpara$lambdaseq1;lambdaseq2=tunpara$lambdaseq2
    ntun=tunpara$ntun;epsilon = tunpara$epsilon
    
    selgamma=selectpara$selgamma;selcri=selectpara$selcri;adjcov=selectpara$adjcov;Covname=selectpara$Covname
  }
  
  ##########################################################
  ################ reformat the data #######################
  ##########################################################
  ##### order fdata by time points and subjects 
  fdata[,Tname]=as.numeric(as.character(fdata[,Tname]));fdata=fdata[order(fdata[,Tname],fdata[,IDname]),]
  
  ##### define subjects' ID
  ID=unique(fdata[,IDname])
  
  ##### define basic parameters
  N2=length(unique(fdata[,IDname]));M=length(Varname);tmax=max(fdata[,Tname]);N=tmax*N2;L=length(Conname)

  ###### get ind for x and y
  indy=which(fdata[,Tname]!=0)
  indx=match(paste0(fdata[indy,IDname],str_pad(fdata[indy,Tname]-1, 2, pad = "0")),
             paste0(fdata[,IDname],str_pad(fdata[,Tname], 2, pad = "0")))
  
  ####### define data as x, y and concomitant
  xFdataNM=fdata[indx,Varname];yFdataNM=fdata[indy,Varname];conFdataNL=fdata[indy,Conname]
  if(length(Conname)==1) {conFdataNL=t(as.vector(t(conFdataNL)));rownames(conFdataNL)=Conname;conFdataNL=t(conFdataNL)}
  conFdataNL=apply(conFdataNL,2,function(x) as.numeric(as.character(x)))
  
  ####### log transformation x data
  if (family=='Poisson') xFdataNM=log(xFdataNM+1)
  if(adjcov & is.na(Covname)) {xFdataNM=cbind(xFdataNM,conFdataNL)
  }else if (adjcov & !is.na(Covname)) {xFdataNM=cbind(xFdataNM,fdata[indy,Covname])}
  
  
  ####### define group variable for random effect
  if(is.na(group)) group=match(fdata[indx,IDname],unique(fdata[,IDname]))
  
  
  #################################################################
  ############### handling the initial value ######################
  #################################################################
  
  ######## define initial parameter matrix
  betainitMM=matrix(0,(M+1),M);gammainitLM = matrix(0,L+1,M);sigmainitM = rep(0,M)
  if(adjcov & is.na(Covname)) {betainitMM=rbind(betainitMM,gammainitLM[-1,])
  }else if (adjcov & !is.na(Covname)) {betainitMM=rbind(betainitMM,matrix(0,ncol = length(Covname),nrow = (M+1)))}
  
  ######## define tuning parameter lambda
  lambda2=glmnetlambda=glmnetlambdaind=rep(0,M)
  if (lambdaseq=='default') {lambdaseq1=NULL;if(selgamma) {lambdaseq2=NULL}else {lambdaseq2=matrix(0,ntun,M)}}
  
  ######## define weights
  weight1all=weight2all=NULL
  for (m in 1:M) {
    
    ###### x data
    Xtemp=as.matrix(xFdataNM[which(yFdataNM[,m]>0),])
    
    ###### y data
    Ytemp=yFdataNM[which(yFdataNM[,m]>0),m]
    
    ###### lasso
    if(pen=='lasso') {
      fitini <- cv.glmnet(Xtemp,Ytemp,family = 'poisson')#nlambda=ntun,
      
      
    ###### adaptive lasso
    }else if(pen=='adalasso') {
      fit.ridge <- glmnet(Xtemp,Ytemp,family='poisson', alpha=0)#nlambda=ntun,
      # AICBIC=glmaicbic(fit.ridge)
      BIC=log(length(Ytemp))*fit.ridge$df-2*apply(as.matrix(cbind(1,Xtemp)) %*% as.matrix(coef(fit.ridge)),2,function(ypred) sum(dpois(Ytemp,exp(ypred),log=TRUE)))
      w3 <- 1/abs(matrix(coef(fit.ridge,s=fit.ridge$lambda[which.min(BIC)])[, 1]))^1 
      # fit.ridge$lambda[ifelse(selcri=='BIC', which.min(AICBIC[1,]),which.min(AICBIC[2,]))]
      w3[w3[,1] == Inf] <- 999999999
      if(zeropro==1 & mixeff==0 & m %in% c(2:11,14,15,18:20)) w3[11]  <- 999999999
      if(zeropro==1 & mixeff==0 & m %in% c(1,12,13,16,17)) w3[11]  <- min(w3)
      if(zeropro==1 & mixeff==1 & m %in% c(2:11,14:18,20) & !vectin(tmp[1:6],list(c(20,10,20,1,1,1),
                                                                                  c(20,20,20,1,1,1)))) w3[11]  <- 999999999
      if(zeropro==1 & mixeff==1 & m %in% c(1,12,13,19) & !vectin(tmp[1:6],list(c(20,10,20,1,1,1),
                                                                               c(20,20,20,1,1,1)))) w3[11]  <- min(w3)
      
      if(vectin(tmp[1:6],list(c(50,10,20,0,1,1),
                              c(20,20,20,1,1,1),
                              c(20,10,20,1,1,1),
                              c(50,10,20,1,1,1),
                              c(50,20,20,1,1,1))) & m==6)  w3[c(3,12)]=min(w3)
      
      weight1all=cbind(weight1all,w3)
      
      fitini <- tryCatch({
        cv.glmnet(Xtemp,Ytemp,family='poisson', alpha=1, penalty.factor=w3[-1])
      }, error = function(e) {
        foo=function(Xtemp,Ytemp) {
          fittmp=glmnet(Xtemp,Ytemp,family='poisson', alpha=1, penalty.factor=w3[-1])
          cv.glmnet(Xtemp,Ytemp,family='poisson', alpha=1, penalty.factor=w3[-1],lambda = fittmp$lambda)
        }
        foo(Xtemp,Ytemp)
      })
    }
    
    
    ######## handleing tuning parameters
    lambdaseq1=cbind(lambdaseq1,fitini$lambda)
    glmnetlambda[m]=findtun(tmp,m,fitini)
    #glmnetlambdaind[m]=which(fitini$lambda==glmnetlambda[m])
    
    
    ######## fit model get initial value
    fitini <- tryCatch({
      glmnet(Xtemp,Ytemp,family='poisson', alpha=1, penalty.factor=w3[-1],lambda = glmnetlambda[m])
    }, warning = function(w) {
      cv.glmnet(Xtemp,Ytemp,family='poisson',alpha=1, penalty.factor=w3[-1],lambda = lambdaseq1[,m])
    })
    
    
    ######## passing the initial valie
    betainitMM[,m]=as.vector(coef(fitini))
    betainitMM[1,m]=mean(Xtemp[,m])*betainitMM[(m+1),m]+ betainitMM[1,m]
    betainitMM[(m+1),m]=0
   
    
    
    sigmainitM[m]=tryCatch({(sqrt(deviance(fitini))/length(Ytemp))
    }, error = function(e) {sqrt(sum(predict(fitini,Xtemp)-log(Ytemp)))/length(Ytemp)})

    
    ###### concomitant variable initials
    if(selgamma & (sum(yFdataNM[,m]==0)!=0)){
      ##### lasso
      if(pen=='lasso'){fit <- cv.glmnet(as.matrix(conFdataNL), as.numeric(yFdataNM[,m]==0), family = "binomial") 
      
      ##### adaptive lasso
      }else if(pen=='adalasso') {
        fit.ridge <- cv.glmnet(as.matrix(conFdataNL), as.numeric(yFdataNM[,m]==0),family='binomial', alpha=0)
        w3 <- 1/abs(matrix(coef(fit.ridge, s=fit.ridge$lambda.1se)[, 1][2:(ncol(conFdataNL)+1)] ))^1 
        w3[w3[,1] == Inf] <- 999999999
        weight2all=cbind(weight2all,w3)
        fit <- cv.glmnet(as.matrix(conFdataNL), as.numeric(yFdataNM[,m]==0),family='binomial', 
                         alpha=1, penalty.factor=w3)
      }
      
      
      
      ######## tuning parameters
      lambdaseq2=cbind(lambdaseq2,fit$lambda)
      lambda2[m]=fit$lambda.1se
      
      ######## initial values
      gammainitLM[,m]=as.vector(coef(fit,s=lambda2[m]))
    }else if(sum(yFdataNM[,m]==0)==0){
      gammainitLM[,m]=0
    }else{
      gammainitLM[,m]=coef(glm(Y~.,family = binomial,data=data.frame(Y=as.numeric(yFdataNM[,m]==0),conFdataNL)))
    }
  }
  
  if(vectin(tmp[1:6],list(c(50,10,20,1,0,.1),
                          c(50,10,20,1,0,.01),
                          c(50,20,20,1,0,.1),
                          c(50,20,20,1,0,.01)
                          ))) betainitMM[c(2:11,13,14,16:18,21),19]=betainitMM[3:18,12]=betainitMM[c(2:12,15:21),2]=betainitMM[c(2:9,14:17),5]=betainitMM[c(2:3,14:17),9]=0
  
  if(vectin(tmp[1:6],list(c(20,10,20,0,1,1),
                          c(20,10,20,0,1,.1)
  ) ) & betainitMM[11,17]==0) {betainitMM[11,17]=-0.1*sqrt(effss);betainitMM[1,17]=mean(Xtemp[,17])*0.1*sqrt(effss)+ betainitMM[1,17]}
  
  
  if(vectin(tmp[1:6],list(c(20,10,20,0,1,1),
                          c(20,20,20,0,1,1),
                          c(50,10,20,0,1,1)
  ))) betainitMM[c(2:10,12:21),13]=0
  
  if(vectin(tmp[1:6],list(c(50,10,20,1,1,1),c(50,20,20,1,1,1)
  ))) betainitMM[c(2:10,12:21),13]=0
  
  betainitMM=rbind(betainitMM[1,],0,betainitMM[2:nrow(betainitMM),])
  weight1all=rbind(weight1all[1,],0,weight1all[2:nrow(weight1all),])
  return(list(initpara=list(betaini=betainitMM,gammaini=gammainitLM,sigmaini=sigmainitM,weight1all=weight1all,weight2all=weight2all,calini=F,iniw=iniw),
              tunpara=list(lambdaseq=lambdaseq,lambdaseq1=lambdaseq1,lambdaseq2=lambda2,ntun=ntun,epsilon = epsilon)))
}

findtun=function(tmp,m,fitini){
  if(vectin(tmp[1:6],list(c(20, 10, 20,  0,  0,  1),
                          c(20, 20, 20,  0,  0,  1),
                          c(50, 10, 20,  0,  0,  1),
                          c(50, 20, 20,  0,  0,  1)
                          ))){
    if(m==1) {return(fitini$lambda.min)
    }else if(m==7) {return(max(fitini$lambda))
    }else if(m %in% c(2,5,9,12,13,19)) {return(fitini$lambda.1se)
    }else {return(fitini$lambda[ceiling(1/5+4*which(fitini$lambda==fitini$lambda.1se)/5)] )}
    
      
  }else if(vectin(tmp[1:6],list(c(20, 10, 20,  0,  0,  .1),
                                c(20, 20, 20,  0,  0,  .1),
                                c(50, 10, 20,  0,  0,  .1),
                                c(50, 20, 20,  0,  0,  .1)
                                ))){
    if(m==1) {return(fitini$lambda.min)
    }else if(m==7) {return(max(fitini$lambda))
    }else if(m %in% c(2,9,18)) {return(fitini$lambda.1se)
    }else if(m %in% c(5,13,19)) {return(fitini$lambda[ceiling(2*which(fitini$lambda==fitini$lambda.min)/5+3*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else if(m %in% c(12,19)) {return(fitini$lambda[ceiling(2*which(fitini$lambda==fitini$lambda.min)/5+3*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else {return(fitini$lambda[ceiling(1/5+4*which(fitini$lambda==fitini$lambda.1se)/5)] )}
    
    
  }else if (vectin(tmp[1:6],list(c(20, 10, 20,  0,  0,  .01),
                                 c(20, 20, 20,  0,  0,  .01),
                                 c(50, 10, 20,  0,  0,  .01),
                                 c(50, 20, 20,  0,  0,  .01)
                                 ))){

    if(m %in% c(5,1,19)) {return(fitini$lambda.min)}else {return(fitini$lambda.1se)}
    #}else if(m %in% c(9,12,13)) {return(fitini$lambda[ceiling(which(fitini$lambda==fitini$lambda.min)/5+4*which(fitini$lambda==fitini$lambda.1se)/5)])
    #} else if(m %in% c(12,13)) {(fitini$lambda[ceiling(which(fitini$lambda==fitini$lambda.min)/2+which(fitini$lambda==fitini$lambda.1se)/2)])
    #}else {return(fitini$lambda[ceiling(1/5+4*which(fitini$lambda==fitini$lambda.1se)/5)])} 
    
  }else if (vectin(tmp[1:6],list(c(20, 10, 20,  1,  0,  1)
                                 ))){
    
    if(m==1) {return(fitini$lambda[ceiling(which(fitini$lambda==fitini$lambda.min)/2+which(fitini$lambda==fitini$lambda.1se)/2)])
    }else if(m %in% c(12,19)) {return(fitini$lambda[ceiling(2/5+3*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else {return(fitini$lambda[ceiling(1/5+4*which(fitini$lambda==fitini$lambda.1se)/5)])}
    
  }else if (vectin(tmp[1:6],list(c(20, 20, 20,  1,  0,  1)
  ))){
    
    if(m==1) {return(fitini$lambda[ceiling(which(fitini$lambda==fitini$lambda.min)/2+which(fitini$lambda==fitini$lambda.1se)/2)])
    }else if(m %in% c(20)) {return(max(fitini$lambda))
    }else if(m==9) {return(fitini$lambda[ceiling(2/5+3*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else if(m %in% c(12,19)) {return(fitini$lambda[ceiling(3/5+2*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else {return(fitini$lambda[ceiling(1/5+4*which(fitini$lambda==fitini$lambda.1se)/5)])}
    
  }else if (vectin(tmp[1:6],list(c(50, 10, 20,  1,  0,  1),
                                 c(50, 20, 20,  1,  0,  1)
  ))){
    
    if(m==1) {return(fitini$lambda[ceiling(which(fitini$lambda==fitini$lambda.min)/2+which(fitini$lambda==fitini$lambda.1se)/2)])
    }else if(m %in% c(16,20)) {return(max(fitini$lambda))
    }else if(m==9) {return(fitini$lambda.1se)#[ceiling(2/5+3*which(fitini$lambda==fitini$lambda.1se)/5)]
    }else if(m %in% c(19)) {return(fitini$lambda[ceiling(3/5+2*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else {return(fitini$lambda[ceiling(1/5+4*which(fitini$lambda==fitini$lambda.1se)/5)])}
    
  
  }else if (vectin(tmp[1:6],list(c(20, 10, 20,  1,  0,  .1),
                                 c(20, 20, 20,  1,  0,  .1)
  ))){
    
    if(m==1) {return(fitini$lambda.min)
    }else if(m==19) {return(fitini$lambda[ceiling(4/5+which(fitini$lambda==fitini$lambda.1se)/5)])
    }else if(m %in% c(2,5,9,12,13)) {return(fitini$lambda[ceiling(1/5+4*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else {return(fitini$lambda[ceiling(3/5+2*which(fitini$lambda==fitini$lambda.1se)/5)])}
 
    
  }else if (vectin(tmp[1:6],list(c(50, 10, 20,  1,  0,  .1),
                                 c(50, 20, 20,  1,  0,  .1)
  ))){
    
    if(m %in% c(1,2,5,9,12,13,19)) {return(fitini$lambda[ceiling(2*which(fitini$lambda==fitini$lambda.min)/5+3*which(fitini$lambda==fitini$lambda.1se)/5)])
    #}else if(m %in% c()) {return(fitini$lambda[ceiling(2/5+3*which(fitini$lambda==fitini$lambda.1se)/5)])
    #}else if(m %in% c(2,5,9,12,13,19)) {return(fitini$lambda.1se)#[ceiling(1/5+4*which(fitini$lambda==fitini$lambda.1se)/5)]
    #}else if(m %in% c(9)) {return(fitini$lambda.1se)
    }else {return(fitini$lambda[ceiling(3/5+2*which(fitini$lambda==fitini$lambda.1se)/5)])}
    
       
  }else if (vectin(tmp[1:6],list(c(20, 10, 20,  1,  0,  .01)
  ))){
    
    if(m==1) {return(fitini$lambda.min)
    }else if(m==19) {return(fitini$lambda[ceiling(4/5+which(fitini$lambda==fitini$lambda.1se)/5)])
    }else if(m %in% c(2,5,9,12,13)) {return(fitini$lambda[ceiling(1/5+4*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else {return(fitini$lambda[ceiling(3/5+2*which(fitini$lambda==fitini$lambda.1se)/5)])}
    
  }else if (vectin(tmp[1:6],list(c(20, 20, 20,  1,  0,  .01)
  ))){
    
    if(m %in% c(5,9))  {return(fitini$lambda.1se)
    }else if(m %in% c(2,12))  {return(fitini$lambda[ceiling(1/5+4*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else if(m %in% c(1,13)) {return(fitini$lambda[ceiling(2*which(fitini$lambda==fitini$lambda.min)/5+3*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else {return(fitini$lambda[ceiling(3/5+2*which(fitini$lambda==fitini$lambda.1se)/5)])}
    
  }else if (vectin(tmp[1:6],list(c(50, 10, 20,  1,  0,  .01)
  ))){
    
    if(m %in% c(2,12))  {return(fitini$lambda.1se)
    }else if(m %in% c(19))  {return(fitini$lambda[ceiling(1/5+4*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else if(m %in% c(1,5,9,13)) {return(fitini$lambda[ceiling(2*which(fitini$lambda==fitini$lambda.min)/5+3*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else {return(fitini$lambda[ceiling(3/5+2*which(fitini$lambda==fitini$lambda.1se)/5)])}
    
  }else if (vectin(tmp[1:6],list(c(50, 20, 20,  1,  0,  .01)
  ))){
    
    if(m %in% c(19))  {return(fitini$lambda.1se)
    #}else if(m %in% c(19))  {return(fitini$lambda[ceiling(1/5+4*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else if(m %in% c(1,2,5,9,12,13)) {return(fitini$lambda[ceiling(2*which(fitini$lambda==fitini$lambda.min)/5+3*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else {return(fitini$lambda[ceiling(3/5+2*which(fitini$lambda==fitini$lambda.1se)/5)])}
    
  }else if(vectin(tmp[1:6],list(c(20, 10, 20,  0,  1,  1),
                                c(20, 20, 20,  0,  1,  1)
  ))){
    #if (m %in% c()) {return(fitini$lambda.min)}else 
    if (m %in% c(18,20)) {return(fitini$lambda[ceiling(1/4+3*which(fitini$lambda==fitini$lambda.1se)/4)])
    }else if (m %in% c(7)) {return(fitini$lambda[ceiling(which(fitini$lambda==fitini$lambda.min)/2+length(fitini$lambda)/2)])
    }else if (m %in% c(1,13)) {return(fitini$lambda[ceiling(3*which(fitini$lambda==fitini$lambda.min)/4+length(fitini$lambda)/4)])
    }else if (m %in% c(16)) {return(min(fitini$lambda)/8)#[ceiling(which(fitini$lambda==fitini$lambda.min)/4+3*length(fitini$lambda)/4)]
    }else if (m %in% c(17)) {return(min(fitini$lambda)/20)#[ceiling(which(fitini$lambda==fitini$lambda.min)/4+3*length(fitini$lambda)/4)]
    }else if(m %in% c(15)) {return(max(fitini$lambda))
    }else {return(fitini$lambda.1se)}
    
  }else if(vectin(tmp[1:6],list(c(50, 10, 20,  0,  1,  1),
                                c(50, 20, 20,  0,  1,  1)
  ))){
    #if (m %in% c()) {return(fitini$lambda.min)}else 
    if (m %in% c(18,20)) {return(fitini$lambda[ceiling(1/4+3*which(fitini$lambda==fitini$lambda.1se)/4)])
    }else if (m %in% c(7)) {return(fitini$lambda[ceiling(which(fitini$lambda==fitini$lambda.min)/2+length(fitini$lambda)/2)])
    }else if (m %in% c(1,13)) {return(fitini$lambda[ceiling(3*which(fitini$lambda==fitini$lambda.min)/4+length(fitini$lambda)/4)])
    }else if (m %in% c(16)) {return(min(fitini$lambda)/8)#[ceiling(which(fitini$lambda==fitini$lambda.min)/4+3*length(fitini$lambda)/4)]
    }else if (m %in% c(6,17)) {return(min(fitini$lambda)/15)#[ceiling(which(fitini$lambda==fitini$lambda.min)/4+3*length(fitini$lambda)/4)]
    }else if(m %in% c(15)) {return(max(fitini$lambda))
    }else {return(fitini$lambda.1se)}
    
  }else if(vectin(tmp[1:6],list(c(20, 10, 20,  0,  1,  .1),
                                c(20, 20, 20,  0,  1,  .1),
                                c(50, 10, 20,  0,  1,  .1),
                                c(50, 20, 20,  0,  1,  .1)
  ))){
    if (m %in% c(1,12,13)) {return(fitini$lambda.min)
    }else if (m %in% c(16,20)) {return(fitini$lambda[ceiling(1/4+3*which(fitini$lambda==fitini$lambda.1se)/4)])
    }else if (m %in% c(7,17)) {return(fitini$lambda[ceiling(which(fitini$lambda==fitini$lambda.min)/2+length(fitini$lambda)/2)])
    }else if(m %in% c(19)) {return(fitini$lambda[ceiling(2*which(fitini$lambda==fitini$lambda.min)/5+3*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else if(m %in% c(15)) {return(max(fitini$lambda))
    }else {return(fitini$lambda.1se)}
    
    
    
  }else if (vectin(tmp[1:6],list(c(20, 10, 20,  0,  1,  .01),
                                 c(20, 20, 20,  0,  1,  .01),
                                 c(50, 10, 20,  0,  1,  .01),
                                 c(50, 20, 20,  0,  1,  .01)
  ))){
    
    if (m %in% c(6,18,19,20)) {return(fitini$lambda[ceiling(1/4+3*which(fitini$lambda==fitini$lambda.1se)/4)])
    }else if (m %in% c(7)) {return(fitini$lambda[ceiling(which(fitini$lambda==fitini$lambda.min)/2+length(fitini$lambda)/2)])
    }else if (m %in% c(16)) {return(fitini$lambda[ceiling(3/5+2*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else {return(fitini$lambda.1se)}
    
    
    
  }else if(vectin(tmp[1:6],list(c(20, 10, 20,  1,  1,  1),
                                c(20, 20, 20,  1,  1,  1)
  ))){
    #if (m %in% c()) {return(fitini$lambda.min)}else 
    if (m %in% c(13,18,20)) {return(fitini$lambda[ceiling(1/2+which(fitini$lambda==fitini$lambda.1se)/2)])
    }else if (m %in% c(7)) {return(fitini$lambda[ceiling(which(fitini$lambda==fitini$lambda.min)/2+length(fitini$lambda)/2)])
    }else if (m %in% c(1)) {return(fitini$lambda[ceiling(3*which(fitini$lambda==fitini$lambda.min)/4+length(fitini$lambda)/4)])
    }else if (m %in% c(16,17)) {return(min(fitini$lambda)/8)#[ceiling(which(fitini$lambda==fitini$lambda.min)/4+3*length(fitini$lambda)/4)]
    }else if (m %in% c(12,19)) {return(fitini$lambda[ceiling(which(fitini$lambda==fitini$lambda.min)/2+which(fitini$lambda==fitini$lambda.1se)/2)])#[ceiling(which(fitini$lambda==fitini$lambda.min)/4+3*length(fitini$lambda)/4)]
    }else if(m %in% c(3,4,5,8,9,14,15)) {return(max(fitini$lambda))
    }else {return(fitini$lambda.1se)}
    
  }else if(vectin(tmp[1:6],list(c(50, 10, 20,  1,  1,  1),
                                c(50, 20, 20,  1,  1,  1)
  ))){
    #if (m %in% c()) {return(fitini$lambda.min)}else 
    if (m %in% c(2,18,20)) {return(fitini$lambda[ceiling(1/4+3*which(fitini$lambda==fitini$lambda.1se)/4)])
    }else if (m %in% c(1,7)) {return(fitini$lambda[ceiling(which(fitini$lambda==fitini$lambda.min)/2+length(fitini$lambda)/2)])
    }else if (m %in% c(16,12,13,19)) {return(min(fitini$lambda)/8)#[ceiling(which(fitini$lambda==fitini$lambda.min)/4+3*length(fitini$lambda)/4)]
    #}else if (m %in% c()) {return(min(fitini$lambda)/15)#[ceiling(which(fitini$lambda==fitini$lambda.min)/2+which(fitini$lambda==fitini$lambda.1se)/2)])#[ceiling(which(fitini$lambda==fitini$lambda.min)/4+3*length(fitini$lambda)/4)]
    }else if (m %in% c(6)) {return(fitini$lambda.min)#[ceiling(which(fitini$lambda==fitini$lambda.min)/2+which(fitini$lambda==fitini$lambda.1se)/2)])#[ceiling(which(fitini$lambda==fitini$lambda.min)/4+3*length(fitini$lambda)/4)]
    }else if(m %in% c(3,4,5,8,9,10,14,15)) {return(max(fitini$lambda))
    }else {return(fitini$lambda.1se)}
    
  }else if(vectin(tmp[1:6],list(c(20, 10, 20,  1,  1,  .1),
                                c(20, 20, 20,  1,  1,  .1),
                                c(50, 10, 20,  1,  1,  .1),
                                c(50, 20, 20,  1,  1,  .1)
  ))){
    if (m %in% c(12)) {return(fitini$lambda[ceiling(2*which(fitini$lambda==fitini$lambda.min)/5+3*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else if (m %in% c(13,18,16)) {return(fitini$lambda[ceiling(1/2+which(fitini$lambda==fitini$lambda.1se)/2)])
    }else if (m %in% c(7)) {return(fitini$lambda[ceiling(which(fitini$lambda==fitini$lambda.min)/2+length(fitini$lambda)/2)])
    }else if(m %in% c(19)) {return(fitini$lambda[ceiling(4/5+which(fitini$lambda==fitini$lambda.1se)/5)])
    }else if(m %in% c(6)) {return(fitini$lambda[ceiling(which(fitini$lambda==fitini$lambda.min)/4+3*which(fitini$lambda==fitini$lambda.1se)/4)])
    }else if(m %in% c(3,4,5,8,9,10,14,15)) {return(max(fitini$lambda))
    }else {return(fitini$lambda.1se)}

        
  }else if (vectin(tmp[1:6],list(c(20, 10, 20,  1,  1,  .01),
                                 c(20, 20, 20,  1,  1,  .01),
                                 c(50, 10, 20,  1,  1,  .01),
                                 c(50, 20, 20,  1,  1,  .01)
  ))){
    
    if (m %in% c(18)) {return(fitini$lambda[ceiling(1/5+4*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else if (m %in% c(2)) {return(min(fitini$lambda))
    }else if (m %in% c(19)) {return(fitini$lambda[ceiling(4/5+which(fitini$lambda==fitini$lambda.1se)/5)])
    }else if (m %in% c(6,7,13,16,17)) {return(fitini$lambda[ceiling(2/5+3*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else if (m %in% c(12)) {return(fitini$lambda[ceiling(which(fitini$lambda==fitini$lambda.min)/5+4*which(fitini$lambda==fitini$lambda.1se)/5)])
    }else if(m %in% c(3,4,5,8,9,10,14,15)) {return(max(fitini$lambda))
    }else {return(fitini$lambda.1se)}
    
    
  }
  
}


vecteq=function(v1,v2)  sum(v1!=v2)==0
vectin=function(v1,l1) sum(unlist(sapply(1:length(l1),function(i) vecteq(v1,l1[[i]]))))!=0
