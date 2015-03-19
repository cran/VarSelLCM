\name{VSLCMpartitions-class}
\docType{class}
\alias{VSLCMpartitions-class}

\title{Class \code{"VSLCMpartitions"}}
\description{
This class contains the estimated partitions
}



\section{Slots}{
  \describe{
    \item{\code{zMAP}:}{A vector indicating the class membership of each individual by using the MAP rule computed for the best model with its maximum likelihood estimates}
    \item{\code{zOPT}:}{Partition maximizing the integrated complete-data likelihood of the selected model}
    \item{\code{tik}:}{Fuzzy partition computed for the best model with its maximum likelihood estimates}
  }
}


\examples{
showClass("VSLCMpartitions")
}

\keyword{classes}
