\name{VSLCMcriteria-class}
\docType{class}
\alias{VSLCMcriteria-class}

\title{Class \code{"VSLCMcriteria"}}
\description{
This class contains the information criteria
}



\section{Slots}{
  \describe{
    \item{\code{loglikelihood}:}{Value of the log-likelihood}
    \item{\code{BIC}:}{Values of the BIC criterion}
    \item{\code{ICL}:}{Values of the ICL criterion}
    \item{\code{MICL}:}{Values of the MICL criterion}  
    \item{\code{degeneracyrate}:}{Rate of EM having degenerated}  
  }
}


\examples{
showClass("VSLCMcriteria")
}

\keyword{classes}
