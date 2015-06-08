\docType{methods}
\name{print}
\alias{print}
\alias{print,VSLCMresultsContinuous-method}
\title{Print an object of class VSLCMresults}
\arguments{
  \item{x}{An object of class \code{\linkS4class{VSLCMresultsContinuous}}}
}
\value{
  NULL. 
}
\description{
 Print an object of class VSLCMresultsContinuous
}
\examples{
\dontrun{
require(VarSelLCM)
data(banknote)
results <- VarSelCluster(banknote[,-1], 2, initModel=40)
print(results)
}
}
\seealso{
  \code{\link{print}}
}

