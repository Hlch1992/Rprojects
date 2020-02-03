\name{parameters}
\alias{parameters}
\docType{data}
\title{
 A list of parameters for simulation
}
\description{
This data is a list object which contains the baseline OTU table, sample data, time number and parameters beta matrix gamma matrix sigma and random effects. Details see reference.
}
\usage{data("parameters")}
\format{
  A list with 8 objects as following.
  \describe{
    \item{\code{baseFdataN2M}}{a numeric matrix containning baseline OTU table}
    \item{\code{conFdataN2L}}{a numeric matrix containning time-independent concomitant variables}
    \item{\code{timeN}}{a numeric indicating total number of time points to simulate}
    \item{\code{interceptM}}{a numeric vector indicating intercept of non-zero autoregressive model}
    \item{\code{betaMM}}{a numeric matrix indicating the interaction matrix}
    \item{\code{gammaLM}}{a numeric matrix indicating fixed effects for the zero state logit model}
    \item{\code{sigmaM}}{a numeric vector indicating the measurement error standard deviation for the non-zero auto-regressive model}
    \item{\code{biN2M}}{a numeric vector indicating the random effects part of the non-zero auto-regressive model}

  }
}
\examples{
data(parameters)
}
\keyword{datasets}
