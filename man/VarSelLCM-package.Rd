\name{VarSelLCM-package}
\alias{VarSelLCM-package}
\alias{VarSelLCM}
\docType{package}


\title{
 Variable Selection in model-based clustering managed by the Latent Class Model for analyzing continuous data with missing values.
}

\description{
The R package VarSelLCM uses a finite mixture model for performing the cluster analysis with variable selection of continuous data by assuming independence between classes. The package deals dataset with missing values by assuming that values are missing at random. The one-dimensional marginals of the components follow Gaussian distributions for facilitating both model interpretation and model selection. The variable selection is led by an alternated optimization procedure for maximizing the MICL criterion. The maximum likelihood inference is done by an EM algorithm for the selected model. This package also performs the imputation of  missing values.
}

\details{
\tabular{ll}{
Package: \tab VarSelLCM\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2015-06-08\cr
License: \tab GPL (>= 2)\cr
}
The main functions to use are \code{\link{VarSelCluster}} and \code{\link{VarSelImputation}}.

Function \code{\link{VarSelCluster}} carries out the model selection by maximizing the MICL criterion, then it performs the maximum likelihood estimation of the selected model via an EM algorithm.

Function \code{\link{VarSelImputation}} performs the imputation of  missing values by taking the expectation of the missing values conditionally on the model, its parameters and on the observed variables.

Tool methods \code{summary} and \code{print} are available for facilitating the interpretation.
}

\author{
Matthieu Marbac and Mohammed Sedki
Maintainer: Mohammed Sedki <mohammed.sedki@u-psud.fr>
}

\references{
  M. Marbac and M. Sedki (2015).
  Variable Selection for Model-Based Clustering using the Integrated Completed-Data Likelihood.
  \emph{Preprint}. 
}


\examples{
\dontrun{
require(VarSelLCM)
data(banknote)

results <- VarSelCluster(banknote[,-1], 2, nbcores=2, initModel=40)

summary(results)

print(results)
}
}
\keyword{ package }