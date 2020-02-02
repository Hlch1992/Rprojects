\name{ARZIMM-package}
\alias{ARZIMM-package}
\alias{ARZIMM}
\docType{package}
\title{
\packageTitle{ARZIMM: Statistical modeling and inference of microbial interaction and stability for microbial
dynamical systems}
}
\description{
\packageDescription{This package is developed to model microbial dynamical systems from longitudinal microbiome data and infer microbial interaction and stability.}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{ARZIMM}
\packageIndices{ARZIMM}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
\packageAuthor{Linchen He}

Maintainer: \packageMaintainer{Linchen He}
}
\references{
Modeling and inference of microbial interaction and stability on the longitudinal microbiome data
}

\keyword{microbial dynamic systems; interaction; stability; longitudinal microbiome data}
\examples{
\dontrun{
require(ARZIMM),
data(sampleparameters)
set.seed(1234)
simMixTime(baseFdataN2M=parameters$baseFdataN2M,conFdataN2L=parameters$conFdataN2L,timeN=parameters$timeN,
interceptM=parameters$interceptM,betaMM=parameters$betaMM,gammaLM=parameters$gammaLM,
sigmaM=parameters$sigmaM,biN2M=parameters$biN2M)

data(sampledata)
Varname=colnames(sampledata)[1:20]
Conname=colnames(sampledata)[21:26]
Tname=colnames(sampledata)[27]
IDname=colnames(sampledata)[28]
ARZIMM::ARZIMM(Varname = Varname,Conname = Conname,fdata = sampledata,IDname = IDname,Tname = Tname,
bootpara=list(bootpval=T,nboot=100))

}
}
