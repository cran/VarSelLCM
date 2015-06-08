\name{VarSelImputation}
\alias{VarSelImputation}
\title{VarSelImputation}

\description{
This function performs the imputation on individuals having missing values.
It uses the parameters resulting from function \link{VarSelCluster}.
Imputation is made by taking the expectation of the missing variables conditionally on 
the observed variables and on the parameters.
}
\usage{
VarSelImputation(obj, ind)
}

\arguments{
  \item{obj}{An object resulting from function \link{VarSelCluster} (i.e. an object of
  VSLCMresultsContinuous).}
  
  \item{ind}{
    An vector of indicating which individuals are imputated (by default, all the individuals having missing values are imputated).
  }
  
  
}
\value{
  A data.frame contaning the imputation made on the selected individuals.
}

\examples{
\dontrun{
require(VarSelLCM)
data(banknote)
banknote[1:3,2:3] <- NA

results <- VarSelCluster(banknote[,-1], 2, initModel=40)


ximput <- VarSelImputation(results, c(1:3))

# Individuals without imputation
print(banknote[c(1:3),-1])

# Individuals with imputation
print(ximput)
}
}