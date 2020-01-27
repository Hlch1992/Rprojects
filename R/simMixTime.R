#' A simulation function
#'
#' Simulation for the ARZIMM model.
#' @param baseFdataN2M a matrix of absolute counts at baseline time
#' @param conFdataN2L a matrix of concomitant variables
#' @param timeN a vector of the number of time points for each subjects
#' @param interceptM a vector of the intercepts of the non-zero autoregression model
#' @param  betaMM a matrix of network table
#' @param  gammaLM a matrix of the zero state logit model
#' @param  sigmaM a vector of the measurement error standard deviation for the non-zero autoregression model
#' @param  biN2M a vector of the random effects for the non-zero autoregression model
#' @param  family a character string indicating the distribtuion. default is Poisson
#' @return a list of simulations:
#'   \item{otu.tab}{a data frame of OTU table}
#'   \item{convar}{a data frame of concomitant variables}
#'   \item{subjid}{a vector of subject-time IDs}
#'   \item{df}{the number of non-zero parameter estimates}
#'   \item{zi}{a list of the random effects}
#' @keywords Simulation
#' @export
#' @examples
#' simMixTime()


simMixTime=function(baseFdataN2M,conFdataN2L,timeN,interceptM,betaMM,gammaLM=NULL,sigmaM,biN2M=NULL,family='Poisson') {

  ######## total number of individual; microbiome
  N2=nrow(baseFdataN2M);M=ncol(baseFdataN2M);L=ncol(conFdataN2L);tmax=max(timeN)

  ######## random intercept for each individual
  if(!is.null(biN2M)) {biN2M=apply(log(baseFdataN2M+1),2, function(x) x-mean(x[which(x!=0)]))
  biN2M[which(biN2M<(-5))]=0}

  ########
  if(is.null(gammaLM)) {
    pN2Mt=rep(list(matrix(0,N2,M)),tmax)
  }else{
    pN2Mt=rep(list(matrix(NA,N2,M)),tmax)
    for (t in 1:tmax)   {
      pN2Mt[[t]]=1/(1+exp(-(cbind(1,conFdataN2L) %*% gammaLM)))#[((t-1)*N2+1):(t*N2),]
      pN2Mt[[t]][,which(colSums(gammaLM)==0)]=0
    }
  }


  ########
  lambdaN2Mt=uN2Mt=rep(list(matrix(NA,N2,M)),tmax)
  xtemp=yFdataNM=as.matrix(baseFdataN2M)
  conFdataNL=cbind(conFdataN2L,time=0)
  for(t in 1:tmax) {
    # xtemp[xtemp==0]=exp(1e-6)
    xtemp=log(xtemp+1)
    # xtemp=cbind(xtemp,conFdataNL[((t-1)*N2+1):(t*N2),])
    lambdaN2Mt[[t]]=exp((t(interceptM) %x% rep(1,N2))+biN2M+(xtemp %*% betaMM))
    uN2Mt[[t]]=matrix(rbinom(N2*M,1,pN2Mt[[t]]),nrow=N2)
    xtemp=(1-uN2Mt[[t]])*matrix(rpois(N2*M,lambdaN2Mt[[t]]),nrow = N2)
    yFdataNM=rbind(yFdataNM,xtemp)
    conFdataNL=rbind(conFdataNL,cbind(conFdataN2L,time=t))
  }
  colnames(yFdataNM)=colnames(baseFdataN2M)
  rownames(yFdataNM)=rownames(conFdataNL)=paste(rep(rownames(baseFdataN2M),(tmax+1)),rep(0:tmax,each=N2),sep = '_')
  list(otu.tab=as.data.frame(yFdataNM),convar=as.data.frame(conFdataNL),subjid=rep(rownames(baseFdataN2M),(tmax+1)),zi=uN2Mt)
}
