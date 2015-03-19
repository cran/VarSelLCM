\name{VarSelLCM-package}
\alias{VarSelLCM-package}
\alias{VarSelLCM}
\docType{package}


\title{
Variable Selection for Model-Based Clustering using the Integrated Complete-Data Likelihood of a Latent Class Model
}

\description{
Variable Selection for model-based clustering by using a mixture model of Gaussian 
distributions assuming conditional independence between variables. The algorithm 
carries out the model selection by optimizing the MICL criterion which has a closed 
form for such a distribution.
}

\details{
\tabular{ll}{
Package: \tab VarSelLCM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-11-25\cr
License: \tab GPL (>= 2)\cr
}
The main functions to use is \code{\link{VarSelCluster}}. It carries out the model selection by maximizing the MICL criterion.
}

\author{
Matthieu Marbac and Mohammed Sedki
Maintainer: Mohammed Sedki <mohammed.sedki@u-psud.fr>
}

\references{
  M. Marbac and M. Sedki (2015).
  Variable selection for model-based clustering using the integrated completed-data likelihood.
  \emph{Preprint}. 
}


\examples{
# Variable selection
require(mclust)
require(VarSelLCM)
data(banknote)

results <- VarSelCluster(banknote[,-1], 2, parallel=FALSE)

summary(results)
print(results)

\dontrun{
## To install multitest package : uncomment the following lines
## source("http://bioconductor.org/biocLite.R")
## biocLite("multtest")
require(mclust)
require(VarSelLCM)
require(multtest)
data(golub)

results <- VarSelCluster(t(golub), 2)

summary(results)

}

}

\keyword{ package }
