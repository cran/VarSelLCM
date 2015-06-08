\name{VarSelCluster}
\alias{VarSelCluster}
\title{VarSelCluster}

\description{
Function VarSelCluster performs the cluster analysis of a continuous data set with missing values.
If option vbleSelec=TRUE, then a variable selection is performed by an iterative algorithm to provide
the model maximizing the MICL criterion.
If option paramEstim=TRUE, then the maximum likelihood inference is performed by an EM algorithm for the
selected model.
Missing data are managed by assuming that data are missing at random.
  }
\usage{
VarSelCluster(x, g, initModel=50, vbleSelec=TRUE, 
paramEstim=TRUE, nbcores=1, nbSmall=250, iterSmall=20,
nbKeep=50, iterKeep=10**3, tolKeep=10**(-3))
}

\arguments{
  \item{x}{A data frame or a matrix where rows correspond to observations and columns correspond to variables. Variables must be  "numeric".}
  
  \item{g}{
    An integer specifying the number of mixture components.
  }
  
  \item{initModel}{
    The number of random initializations of the algorithm. 
    The default is 50. 
  }
  
  \item{vbleSelec}{
    Logical indicating if a variable selection is performed. 
    The default is TRUE. 
  }

  
  \item{paramEstim}{
    Logical indicating if the maximum likelihood (ML) inference is performed for the best model (according to the MICL criterion). 
    The default is TRUE. 
  }
  
  \item{nbcores}{
    Number of cores.
    The default is 1.
  }
  \item{nbSmall}{
    Numeric indicating the number of SmallEM algorithms performed for the ML inference. 
    The default is 250. 
  }
  \item{iterSmall}{
    Numeric indicating the number of iterations for each SmallEM algorithm.
    The default is 20. 
  }
  \item{nbKeep}{
    Numeric indicating the number of chains used for the final EM algorithm.
    The default is 50. 
  }
  \item{iterKeep}{
    Numeric indicating the maximal number of iterations for each EM algorithm.
    The default is 1000. 
  }
  \item{tolKeep}{
    Numeric indicating the maximal gap between two successive iterations of EM algorithm which stops the algorithm.
    The default is 0.001. 
  }

}
\value{
  An object of class \code{\linkS4class{VSLCMresultsContinuous}}  containing the model maximizing the MICL criterion.
  Thus, it indicates which variables are relevant to the clustering. It also provides the maximum likelihood estimates associated to this model.
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

# Cluster analysis without variable selection 
resultswithout <- VarSelCluster(banknote[,-1], 2, vbleSelec=FALSE)
summary(resultswithout)

# Cluster analysis with variable selection 
results <- VarSelCluster(banknote[,-1], 2, nbcores=2, initModel=40)
summary(results)
}
}