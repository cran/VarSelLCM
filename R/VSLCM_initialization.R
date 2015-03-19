VSLCM_initialization_omega <- function(n, d, g){
  omega <- rep(0, d)
  cb <- sample(2:(min(n,d)),1)
  omega[sample(1:d, cb)] <- 1
  return(omega)
}



VSLCM_initialization_priors <- function(x){
  priors <- matrix(1, ncol(x), 4)
  colnames(priors) <- c("alpha", "beta", "lambda", "delta")
  for (j in 1:ncol(x)){
    priors[j,1] <- 1.28*2
    priors[j,2] <- sqrt(0.72 * var(x[,j]))
    priors[j,3] <- mean(x[,j])
    priors[j,4] <- 2.6 /(max(x[,j]) - min(x[,j]))
  }
  #priors <- matrix(1, ncol(x), 4)
  return(priors)
}

VarSelStartingPoint <- function(x, g, omega, z, priors){
  if (missing(priors))
    priors <- VSLCM_initialization_priors(as.matrix(x))

  # Initialization
  if (missing(omega))
    omega <- VSLCM_initialization_omega(nrow(x), ncol(x), g)
  
  if (missing(z))
    z <- kmeans(as.matrix(scale(x[,which(omega==1)],TRUE,TRUE)), g)$cluster
  
  if (min(z)>0)
    z <- z-1
  
  if (is.null(colnames(x)))
    colnames(x) <- paste("X",1:ncol(x), sep="")
  
  
  starting <- new("VSLCMresults",
                  data = as.matrix(x),
                  priors = priors,
                  criteria = new("VSLCMcriteria", likelihood=-Inf, BIC=-Inf, ICL=-Inf, MICL=-Inf),
                  partitions = new("VSLCMpartitions", zMAP=z, zOPT=z),
                  model = new("VSLCMmodel", g=g, omega=omega),
                  parameters = new("VSLCMparameters")
  )

  return(starting)
}
