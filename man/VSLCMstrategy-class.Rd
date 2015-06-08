\name{VSLCMstrategy-class}
\docType{class}
\alias{VSLCMstrategy-class}

\title{Class \code{"VSLCMstrategy"}}
\description{
This class contains the estimated partitions
}



\section{Slots}{
  \describe{
    \item{\code{initModel}:}{number of initialisations for the model selection algorithm}
    \item{\code{vbleSelec}:}{logical indicating if the variable selection is performed}
    \item{\code{paramEstim}:}{logical indicating if the parameter estimation is performed}
    \item{\code{parallel}:}{logical indicating if a parallelisation is done}
    \item{\code{nbSmall}:}{number of small EM}
    \item{\code{iterSmall}:}{number of iteration for the small EM}
    \item{\code{nbKeep}:}{number of chains kept for the EM}
    \item{\code{iterKeep}:}{maximum number of iterations for the EM}
    \item{\code{tolKeep}:}{value of the difference between successive iterations of EM stopping the EM}
  }
}


\examples{
showClass("VSLCMstrategy")
}

\keyword{classes}
