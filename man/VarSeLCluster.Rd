\name{VarSelCluster}
\alias{VarSelCluster}
\title{
VarSelCluster
}

\description{
   This function performs the variable selection for a model-based clustering. 
   The data distribution follow a mixture of Gaussian distributions assuming conditional independence between variables. 
   The model selection is carried out by maximizing the MICL criterion.
}
\usage{
VarSelCluster(x, g, nbinit=50, parallel = TRUE)
}

\arguments{
  \item{x}{
    A matrix, or data frame where rows correspond to observations and columns correspond to variables. 

  }
  \item{g}{
    An integer specifying the number of mixture components
  }
  \item{nbinit}{
    The number of random initializations of the algorithm. 
    The default is 50. 
  }
  \item{parallel}{
    A Boolean indicating if the code is parallelized or not.
    The default is TRUE.
  }
}
\value{
  An object of class \code{\linkS4class{VSLCMresults}} providing the model maximizing the MICL criterion.
  Thus, it indicates which variables are relevant to the clustering. It also provides the maximum
  likelihood estimates associated to this model.
}

\references{
  M. Marbac and M. Sedki (2015).
  Variable selection for model-based clustering using the integrated completed-data likelihood.
  \emph{Preprint}. 
}


\examples{
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
