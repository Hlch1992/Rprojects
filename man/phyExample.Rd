\name{phyExample}
\alias{phyExample}
\docType{data}
\title{
A phyloseq-class experiment-level object example data.
}
\description{
This data is a phyloseq-class object which contains the OTU table, sample data, taxonomy table and phylogenetic tree information. Details about the phyloseq class please see the R package "phyloseq".
}
\usage{data("phyExample")}
\format{
  A phyloseq-class experiment-level object
\describe{
\item{\code{otu_table()}}{   OTU Table:         [ 20 taxa and 1050 samples ]

  \code{..@ .Data}: a matrix where each row represents the suject,each column represents the OTU.

  \code{..@ taxa_are_rows}: logi FALSE
}

\item{\code{sample_data()}}{Sample Data:       [ 1050 samples by 8 sample variables ]

A data frame with 1050 observations for six covariates w1...w6 and two variables time and subjid.

\code{w1...w6}: numeric vectors as binary covariates;

\code{time}: a numeric vector indicating the time of sampling;

\code{subjid}: a string vector indicating which subject the sample belonging to.
}
}
}
\examples{
data(phyExample)
}
\keyword{datasets}
