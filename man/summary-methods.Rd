\docType{methods}
\name{summary}
\alias{summary}
\alias{summary,VSLCMresultsContinuous-method}
\title{Produce result summary of a VSLCMresultsContinuous class}
\arguments{
  \item{object}{An object of class \code{\linkS4class{VSLCMresultsContinuous}}}
}
\value{
  NULL. Summaries to standard out.
}
\description{
  Produce result summary of a VSLCMresultsContinuous class
}
\examples{
\dontrun{
require(VarSelLCM)
data(banknote)
results <- VarSelCluster(banknote[,-1], 2, initModel=40)
summary(results)
}
}
\seealso{
  \code{\link{summary}}
}

