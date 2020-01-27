#' EM algorithm
#'
#' EM algorithm used in the ARZIMM model.
#' @param data a list of data with componenets:
#' \describe{
#'   \item{yFdata}{a vector of}
#'   \item{xFdata}{a matrix of variables corresponding to the non-zero auto-regressive model}
#'   \item{conFdata}{a matrix of covariates corresponding to the zero state logit model}
#'   \item{group}{a vector of numbers as group indicator}
#' }
#' @param para a list of parameter estimates:
#' \describe{
#'   \item{beta}{initial value for beta}
#'   \item{gamma}{initial value for gamma}
#'   \item{sigma}{initial value for sigma}
#' }
#' @param weight a vector of observation weightsfor both the non-zero auto-regressive model and
#' the zero state logit model
#' @param family a character string indicating the distribtuion. default is Poisson
#' @param  selgamma logical; should concomitant variables in the zero state logit model be selected
#' @return a list of fits
#'   \item{para}{a list of parameter estimates:
#'   \describe{
#'     \item{beta}{beta estimates}
#'     \item{gamma}{gamma estimates}
#'     \item{sigma}{sigma estimates}
#'   }}
#'   \item{ciestm}{the estimated random effects}
#'   \item{conv}{logical; did the algorithm converged}
#'   \item{df}{the number of non-zero parameter estimates}
#'   \item{bic}{a vector of BIC, AIC, and log likelihood}
#'   \item{mse}{a vector of square root of mean (pearson) standard error}
#'   \item{lambda}{a vector of \code{lambda} sequence}
#' @keywords EM algorithm
#' @export
#' @examples
#' emnr()

emnr=function(data,para,weight,family,selgamma=F){

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

    ###### M-step:
    maximres=Update_p_M(data=data,para=list(beta=beta,gamma=gamma,sigma=sigma),
                        expectres=expectres,weight=weight,selgamma=selgamma)
    lnew=maximres$logl
    # print(lnew)

    beta=maximres$beta;gamma=maximres$gamma;sigma=maximres$sigma
    if(sum(is.na(c(beta,gamma,sigma)))>0 | is.na(lnew) ) {converge=2;break}

    ####### test convergency
    eps=sum(abs(c(beta- betaold,gamma- gammaold,sigma- sigmaold)))#
    if(its2==1) epsold=eps
    if (eps<SMALL | lnew < lold) {converge=1;break
    }else if(eps>5*BIG | eps>epsold*10 | max(abs(sigma))>BIG | min(sigma) < 0 | abs(beta[1]) > 50 | sum(beta)==0 | max(abs(cihatN2))> 50 ) {
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
  maximres=Update_p_M(data=data,para=list(beta=beta,gamma=gamma,sigma=sigma),
                      expectres=expectres,weight,selgamma=selgamma,sel=T,fast = F)
  beta=maximres$beta;gamma=maximres$gamma;sigma=maximres$sigma;logl=maximres$logl
  if(sum(is.na(c(beta,gamma,sigma)))>0)  {beta = betaold;gamma = gammaold;sigma = sigmaold; cihatN2=cihatN2old}

  para=list(beta=beta,gamma=gamma,sigma=sigma,cihatN2=cihatN2)
  ciestm=list(cihatN2=cihatN2,civarN2=civarN2)
  mf=modfit(data,para)
  rmse=mf$rmse;psmse=mf$psmse;ypred=mf$ypred
  aic=-2*sum(logl)+sum(unlist(para)!=0)*log(N)
  bic=-2*sum(logl)+sum(unlist(para)!=0)*log(N)
  if(is.na(aic)) aic=Inf
  if(is.na(bic)) bic=Inf
  return(list(para=para,ciestm=ciestm,conv=converge,df=sum(unlist(para)!=0),bic=c(aic=aic,bic=bic,logl=logl),
              mse=c(rmse,psmse),ypred=ypred,lambda=maximres$lambda))
}
