\name{sampledata}
\alias{sampledata}
\docType{data}
\title{
a data frame containing all variables to be analysized}
\description{
This is a data frame object which contains OTU varalbes, concomitant variables, time variable and subject ID.
}
\usage{data("sampledata")}
\format{
  A data frame with 1050 observations on the following 28 variables.
  \describe{
    \item{\code{M1...M20}}{numeric vectors indicating the absolute abundance of taxa}
    \item{\code{w1...W6}}{numeric vectors indicating the value of concomitant variables}
    \item{\code{time}}{a numeric vector indicating the time of sampling}
    \item{\code{subjid}}{a factor indicating the subject IDs}
  }
}
\examples{
data(sampledata)
}
\keyword{datasets}
