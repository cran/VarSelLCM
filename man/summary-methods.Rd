\docType{methods}
\name{summary}
\alias{summary}
\alias{summary,VSLCMresults-method}
\title{Produce result summary of a VSLCMresults class}
\arguments{
  \item{object}{An object of class \code{\linkS4class{VSLCMresults}}}
}
\value{
  NULL. Summaries to standard out.
}
\description{
  Produce result summary of a VSLCMresults class
}
\examples{
require(mclust)
require(VarSelLCM)
data(banknote)
results <- VarSelCluster(banknote[,-1], 2, parallel=FALSE)
summary(results)
}
\seealso{
  \code{\link{summary}}
}

