#' Auto-Regressive Zero-Inflated Mixed Model
#'
#' This function allows you to fit the ARZIMM.
#' @param Varname a vector of character string indicating the taxa names in the non-zero auto-regressive model
#' @param Conname a vector of character string indicating the concomitant variable names in the zero state logit model
#' @param Covname a vector of character string indicating the covariate names in the non-zero auto-regressive model
#' @param fdata a data frame containing all variables to be analysized
#' @param IDname a character string indicating the subject ID. Default is ID
#' @param Tname a character string indicating the time variable. Default is Time
#' @param group a vector of numbers as group indicator. To be calculated automatically When absent
#' @param family a character string indicating the distribtuion. default is Poisson
#' @param initpara a named list of user-specified initial values:
#' \describe{
#'   \item{betaini}{the matrix of fixed effects for the non-zero auto-regressive model}
#'   \item{gammaini}{the matrix of fixed effects for the zero state logit model}
#'   \item{sigmaini}{the measurement error standard deviation for both the non-zero auto-regressive model and
#'   the zero state logit model}
#'   \item{weight1all}{observation weights for the non-zero auto-regressive model. Default is 1 for each observation}
#'   \item{weight2all}{observation weights for the zero state logit model. Default is 1 for each observation}
#'   \item{calini=T}{logical; should initial parameters be calculated. Default is Ture}
#'   \item{iniw}{logical; should observation weights be calculated according to initial parameters. Default is False}
#'   When this list of initial values does not contain some of these components or contains components
#'   not of the appropriate length, then the default initial values are used instead.
#' }
#' @param tunpara a list of control values with components:
#' \describe{
#'   \item{lambdaseq1}{a user supplied \code{lambda} sequence for the non-zero auto-regressive model.
#'   Typical usage is to have the program compute its own \code{lambda} sequence. Supplying a value of
#'   \code{lambda} overrides this. WARNING: use with care. Avoid supplying a single value for \code{lambda}}
#'   \item{lambdaseq2}{a user supplied \code{lambda} sequence for the  zero state logit model.
#'   Typical usage is to have the program compute its own \code{lambda} sequence. Supplying a value of
#'   \code{lambda} overrides this. WARNING: use with care. Avoid supplying a single value for \code{lambda}}
#'   \item{ntun}{the number of \code{lambda} values - default is 50}
#'   \item{epsilon}{the range of lambda values; default= c(5e4,5)}
#' }
#' @param selectpara a list of control values with components:
#' \describe{
#'   \item{selgamma}{logical; should concomitant variables in the zero state logit model be selected}
#'   \item{selcri}{method to be use for variable selection. Currently three options. The default is BIC.
#'   Other choices are AIC and CV}
#'   \item{pen}{penalty; defualt is 'adalasso'}
#' }
#' @param bootpara a list of control values with components:
#' \describe{
#'   \item{bootpval}{logical; should p value be calculated using bootstrap}
#'   \item{nboot}{the number of bootstrap simulations; default is 500}
#' }
#'
#' @return an object of class \code{"ARZIMMObject"} is returned, which is a list with the ingredients of fit.
#'   \item{nwtable}{the matrix of network table of fit}
#'   \item{mseest}{a list of mean square error:
#'   \describe{
#'     \item{rmseest}{the root square of mean standard error}
#'     \item{rmseest}{the root square of mean pearson standardized error}
#'   }}
#'   \item{paralist}{a list of parameter estimates:
#'   \describe{
#'     \item{beta}{the matrix of fixed effects for the non-zero auto-regressive model}
#'     \item{gamma}{the matrix of fixed effects for the zero state logit model}
#'     \item{sigma}{the measurement error standard deviation for both the non-zero auto-regressive model and
#'     the zero state logit model}
#'     \item{ciest}{the estimated random effects part of both the non-zero auto-regressive model and
#'     the zero state logit model}
#'   }}
#'   \item{runtime}{running time of the program}
#'   \item{datalist}{an object of class \code{"ARZIMMData"}}
#'   \item{resultall}{a list of parameter estimates of the fits with ingredients of the lambda}
#'   \item{bootparapval}{a list of p values obtain via bootstrap with componenets:
#'   \describe{
#'     \item{betapval}{a vector of p values of fixed effects for the non-zero auto-regressive model}
#'     \item{gammapval}{a vector of p values of fixed effects for the zero state logit model}
#'     \item{sigmapval}{a vector of p values of the measurement error standard deviation for both the
#'     non-zero auto-regressive model and the zero state logit model}
#'   }}
#'   \item{tunlist}{the values of parameters used in the fits with components:
#'   \describe{
#'     \item{lambdabeta}{the values of \code{lambda} used in the non-zero auto-regressive model}
#'     \item{lambdagamma}{the values of \code{lambda} used in the zero state logit model}
#'     \item{weight1all}{observation weights used in the non-zero auto-regressive model}
#'   }}
#'   \item{parasetup}{a list of parameters used to initial the ARZIMM program with components:
#'   \describe{
#'     \item{initpara}{a list of initial parameter inputs; if the inputs are absent, default values are included}
#'     \item{lambda}{the values of \code{lambda} used in the fits.}
#'     \item{tunpara}{a list of tunning parameter inputs; if the inputs are absent, default values are included}
#'     \item{selectpara}{a list of selection parameter inputs; if the inputs are absent, default values are included}
#'   }}
#' @keywords ARZIMM
#' @export
#' @examples
#'
#' data(sampledata)
#' ARZIMM::ARZIMM(Varname = Varname,Conname = Conname,fdata = fdata,IDname = IDname,Tname = Tname,
#' bootpara=list(bootpval=T,nboot=100))



ARZIMM=function(Varname=NA,Conname=NA,Covname=NA,fdata=NA,IDname='ID',Tname='Time',group=NA,family='Poisson',
                   initpara=list(betaini=NULL,gammaini=NULL,sigmaini=NULL,weight1all=NULL,weight2all=NULL,calini=T,iniw=F),
                   tunpara=list(lambdaseq1=NULL,lambdaseq2=NULL,ntun=50,epsilon = c(5e4,5)),
                   selectpara=list(selgamma=T,selcri='BIC',pen='adalasso'),
                   bootpara=list(bootpval=T,nboot=500))
{

  # options(warn=-1)
  start_time <- Sys.time()

  cat("Program is running..be patient...")
  ##########################################################
  ################ reformat the data #######################
  ##########################################################

  if(TRUE){
    ##### order fdata by time points and subjects
    fdata[,Tname]=as.numeric(as.character(fdata[,Tname]));fdata=fdata[order(fdata[,IDname],fdata[,Tname]),]

    ##### define subjects' ID
    ID=unique(fdata[,IDname])

    ##### define basic parameters
    N2=length(unique(fdata[,IDname]));M=length(Varname);N=nrow(fdata)
    if(!is.na(Conname[1])) {L=length(Conname)}else{L=0}
    if(!is.na(Covname[1])) {q=length(Covname)}else{q=0}
    #;tmax=max(fdata[,Tname])

    ###### get ind for x and y
    indy=which(fdata[,Tname]!=0)
    indx=match(paste0(fdata[indy,IDname],str_pad(fdata[indy,Tname]-1, 2, pad = "0")),
               paste0(fdata[,IDname],str_pad(fdata[,Tname], 2, pad = "0")))
    ind0=match(fdata[indy,IDname],fdata[which(fdata[,Tname]==0),IDname])
    xFdataNM0=fdata[ind0,Varname]

    ####### define data as y, x, covariates and concomitant
    yFdataNM=fdata[indy,Varname];  xFdataNM=fdata[indx,Varname]
    if(is.na(Conname[1])){conFdataNL=NULL}else{conFdataNL=t(t(fdata[indx,Conname]))}
    if(is.na(Covname[1])){covFdataNq=NULL}else{covFdataNq=t(t(fdata[indx,Covname]))}
    conFdataNL=apply(conFdataNL,2,function(x) as.numeric(as.character(x)))
    if(!is.na(Covname)) covFdataNq=apply(covFdataNq,2,function(x) as.numeric(as.character(x)))

    ####### log transformation x data
    if (family=='Poisson') xFdataNM=log(xFdataNM+1)
    if(!is.na(Covname[1])) {xFdataNMq=cbind(xFdataNM,covFdataNq)}else{xFdataNMq=xFdataNM}

    ####### define group variable for random effect
    if(is.na(group[1])) group=match(fdata[indx,IDname],unique(fdata[,IDname]))

    datalist=list(yFdataNM=yFdataNM,xFdataNM=xFdataNMq,conFdataNL=conFdataNL,group=group,xFdataNM0=xFdataNM0)

  }


  ##########################################################
  ################ set initial value #######################
  ##########################################################

  if(initpara$calini)  {
    iniest=ARZIMM.initcal(datalist=datalist,family=family,iniw=initpara$iniw, selectpara=selectpara)
    initpara=iniest$initpara;tunpara=iniest$tunpara
  }
  betainitMM=initpara$betaini;gammainitLM=initpara$gammaini;sigmainitM=initpara$sigmaini
  weight1all=initpara$weight1all;weight2all=initpara$weight2all;iniw=initpara$iniw

  lambdaseq1=tunpara$lambdaseq1;lambdaseq2=tunpara$lambdaseq2$lambda2;ntun=tunpara$ntun;epsilon=tunpara$epsilon

  selgamma=selectpara$selgamma;selcri=selectpara$selcri;pen=selectpara$pen

  bootpval=bootpara$bootpval;nboot=bootpara$nboot

  ################################################
  ################# penalty ######################
  ################################################
  if(iniw) {  if(pen=='lasso') {weight1all=matrix(1,nrow = nrow(betainitMM),ncol = ncol(betainitMM))
  weight2all=matrix(1,nrow = nrow(gammainitLM),ncol = ncol(gammainitLM))
  } else if (pen=='adalasso') {weight1all=1/abs(betainitMM);weight2all=1/abs(gammainitLM)}
    weight1all[weight1all == Inf] = weight2all[weight2all == Inf] = max(weight1all[weight1all != Inf])*1000
    weight1all=weight1all[-1,];weight2all=weight2all[-1,]
    weight1all[1,]=0
  }


  ################################################
  ################## progress bar ################
  ################################################
  result=c(NULL)

  #start_time <- Sys.time()

  for (m in 1:M) {
    datam=list(yFdata=yFdataNM[,m],xFdata=cbind(T0=xFdataNM0[,m],xFdataNM),conFdata=conFdataNL,group=group)
    param=list(betainit=betainitMM[,m],gammainit=gammainitLM[,m],sigmainit=sigmainitM[c(m,(m+M))])
    weightm=list(w1=weight1all[,m],w2 = weight2all[,m])
    # betaMntun=gammaLntun=sigmantun=bihatN2ntun=bivarN2ntun=convntun=NULL
    restemp = emnr(data = datam, para = param,family=family, weight = weightm,selgamma=selgamma)#,Lm=Lm

    result=c(result,list(restemp))
  }

  nwtable=beta=gamma=ciest=NULL;sigma=rep(0,2*M)

  for (m in 1:M) {
    nwtable=cbind(nwtable,result[[m]]$para$beta[3:(2+M)])
    beta=cbind(beta,result[[m]]$para$beta)
    gamma=cbind(gamma,result[[m]]$para$gamma)
    sigma[c(m,m+M)]=result[[m]]$para$sigma
    ciest=cbind(ciest,result[[m]]$ciestm$cihatN2)
  }
  rownames(nwtable)=colnames(nwtable)=colnames(beta)=names(sigma)[1:M]=names(sigma)[(M+1):(2*M)]=colnames(ciest)=Varname

  rmseest=psmseest=NULL
  for (m in 1:M) {
    rmseest=c(rmseest,result[[m]]$mse[1])
    psmseest=c(psmseest,result[[m]]$mse[2])
  }

  tunlist=list(lambdabeta=unlist(sapply(1:M,function(x) result[[x]]$lambda[1])),
               lambdagamma=unlist(sapply(1:M,function(x) result[[x]]$lambda[2])),
               weight1all=weight1all,weight2all=weight2all,selgamma=selgamma)
  paralist=list(beta=beta,gamma=gamma,sigma=sigma,ciest=ciest)
  if(bootpval){
    cat("Bootstrap is running..be patient...")

    parabootlist=paralist;tunbootlist=tunlist

    # create progress bar
    pb <- txtProgressBar(min = 0, max = nboot, style = 3)

    bootres=NULL
    for (boottime in 1:nboot){
      set.seed(boottime)
      yFdataNMboot=bootsim(datalist=datalist,paralist=parabootlist)
      datalisttmp=datalist;datalisttmp$yFdataNM=yFdataNMboot

      boottmp=bootest(datalist=datalisttmp,paralist=parabootlist,tunlist=tunbootlist)

      bootres=rbind(bootres,boottmp)
      setTxtProgressBar(pb, boottime)

    }
    ################## progress bar
    close(pb)

    bootsd=apply(bootres, 2, function(x) sd(x))
    pointest=c(as.vector(beta),as.vector(gamma),sigma)
    chistat=unlist(sapply(1:length(pointest), function(x) ifelse(bootsd[x]==0,0,(pointest[x]/bootsd[x])^2)))
    bootwtpval=pchisq(chistat,1,lower.tail = F)
    bootparapval=list(betapval=matrix(bootwtpval[1:length(as.vector(beta))],ncol = M),
                      gammapval=matrix(bootwtpval[(1+length(as.vector(beta))):(length(as.vector(beta))+length(as.vector(gamma)))],ncol = M),
                      sigmapval=bootwtpval[(1+length(as.vector(beta))+length(as.vector(gamma))):length(bootwtpval)])
  }else{bootparapval=NULL}


  end_time <- Sys.time()

  print(end_time - start_time)
  ################## return
  return(list(nwtable=nwtable,mseest=list(rmseest,psmseest),paralist=list(beta=beta,gamma=gamma,sigma=sigma,ciest=ciest),
              runtime=end_time - start_time,datalist=datalist,resultall=result,bootparapval=bootparapval,tunlist=tunlist,
              parasetup=list(initpara=iniest$initpara,lambda=unlist(sapply(1:M,function(x) result[[x]]$lambda[1])),
                             tunpara=tunpara,selectpara=selectpara)))
}

