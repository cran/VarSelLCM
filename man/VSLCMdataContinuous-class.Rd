\name{VSLCMdataContinuous-class}
\docType{class}
\alias{VSLCMdataContinuous-class}

\title{Class \code{"VSLCMdataContinuous"}}
\description{
This class contains the data information
}

\section{Slots}{
  \describe{
    \item{\code{n}:}{Number of individuals.}
    \item{\code{d}:}{Number of continuous variables.}
    \item{\code{data}:}{Whole data set.}
    \item{\code{notNA}:}{Logical indicating if the realization is missing.}
    \item{\code{priors}:}{Values of the hyper-parameters. One row corresponds to one variable, and hyperparameters are set as follows: alpha=beta=1, lambda=mean(x[,j]), delta=1/(100g).}
  }
}


\examples{
showClass("VSLCMdataContinuous")
}

\keyword{classes}
