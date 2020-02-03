#' internal ARZIMM parameters
#'
#' View and/or change the factory default parameters in ARZIMM
#'
#' If called with no arguments, \code{ARZIMM.control()} returns a list with the
#' current settings of these parameters. Any arguments included in the call
#' sets those parameters to the new values, and then silently returns. The
#' values set are persistent for the duration of the R session.
#'
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
#'
#' @return A list with named elements as in the argument list
#' @author Linchen He \cr Maintainer: Linchen He
#' \email{Linchen.He@nyulangone.org}
#' @seealso \code{ARZIMM}
#' @keywords model controls
#' @examples
#'
#' ARZIMM.control(bootpara=list(bootpval=FALSE))
#' ARZIMM.control(selectpara=list(selgamma=FALSE,pen='lasso'))
#'
#' @export ARZIMM.control
#'
ARZIMM.control <-
  function (initpara=list(betaini=NULL,gammaini=NULL,sigmaini=NULL,weight1all=NULL,weight2all=NULL,
                          calini=TRUE,iniw=FALSE),
            tunpara=list(lambdaseq1=NULL,lambdaseq2=NULL,ntun=50,epsilon = c(5e4,5)),
            selectpara=list(selgamma=TRUE,selcri='BIC',pen='adalasso'),
            bootpara=list(bootpval=TRUE,nboot=500))  {

    initpara.default=list(betaini=NULL,gammaini=NULL,sigmaini=NULL,weight1all=NULL,weight2all=NULL,
                  calini=TRUE,iniw=FALSE)
    tunpara.default=list(lambdaseq1=NULL,lambdaseq2=NULL,ntun=50,epsilon = c(5e4,5))
    selectpara.default=list(selgamma=TRUE,selcri='BIC',pen='adalasso')
    bootpara.default=list(bootpval=TRUE,nboot=500)

    if(!is.null(initpara$betaini)) initpara.default$betaini=initpara$betaini
    if(!is.null(initpara$gammaini)) initpara.default$gammaini=initpara$gammaini
    if(!is.null(initpara$sigmaini)) initpara.default$sigmaini=initpara$sigmaini
    if(!is.null(initpara$weight1all)) initpara.default$weight1all=initpara$weight1all
    if(!is.null(initpara$weight2all)) initpara.default$weight2all=initpara$weight2all
    if(!is.null(initpara$calini)) initpara.default$calini=initpara$calini
    if(!is.null(initpara$iniw)) initpara.default$iniw=initpara$iniw


    if(!is.null(tunpara$lambdaseq1)) tunpara.default$lambdaseq1=tunpara$lambdaseq1
    if(!is.null(tunpara$lambdaseq2)) tunpara.default$lambdaseq2=tunpara$lambdaseq2
    if(!is.null(tunpara$ntun)) tunpara.default$ntun=tunpara$ntun
    if(!is.null(tunpara$epsilon)) tunpara.default$epsilon=tunpara$epsilon


    if(!is.null(selectpara$selgamma)) selectpara.default$selgamma=selectpara$selgamma
    if(!is.null(selectpara$selgamma)) selectpara.default$selgamma=selectpara$selgamma
    if(!is.null(selectpara$selgamma)) selectpara.default$selgamma=selectpara$selgamma


    if(!is.null(bootpara$bootpval)) bootpara.default$bootpval=bootpara$bootpval
    if(!is.null(bootpara$nboot)) bootpara.default$nboot=bootpara$nboot

    return(list(initpara=initpara.default, tunpara=tunpara.default,
                selectpara=selectpara.default, bootpara=bootpara.default))
  }
