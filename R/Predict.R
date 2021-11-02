proba.post <- function(object, newdata){

  logprob <- matrix(object@param@pi, nrow(newdata), object@model@g, byrow=TRUE)
  for (nom in colnames(newdata)){
    xnotna <- newdata[,which(colnames(newdata)==nom)]
    where <- which(!is.na(xnotna))
    xnotna <- xnotna[where]
    if (nom %in% rownames(object@param@paramContinuous@mu)){
      who <- which(nom == rownames(object@param@paramContinuous@mu))
      for (k in 1:object@model@g) logprob[where,k] <- logprob[where,k] + dnorm(xnotna, object@param@paramContinuous@mu[who,k], object@param@paramContinuous@sd[who,k], log=TRUE)
    }else if (nom %in% rownames(object@param@paramInteger@lambda)){
      who <- which(nom == rownames(object@param@paramInteger@lambda))
      for (k in 1:object@model@g) logprob[where,k] <- logprob[where,k] + dpois(xnotna, object@param@paramInteger@lambda[who,k], log=TRUE)
    }else if (nom %in% names(object@param@paramCategorical@alpha)) {
      who <- which(nom ==  names(object@param@paramCategorical@alpha))
      for (k in 1:object@model@g){
        for (h in 1:ncol(object@param@paramCategorical@alpha[[who]]))
          logprob[where,k] <- logprob[where,k] + log(object@param@paramCategorical@alpha[[who]][k,h] ** (xnotna == colnames(object@param@paramCategorical@alpha[[who]])[h]))
      }}
  }
  prob <- exp(logprob - apply(logprob, 1, max))
  prob/rowSums(prob)
}

setGeneric ( name= "predict",  def = function(object, newdata, type="probability"){ standardGeneric("predict")})
########################################################################################################################
## predict
########################################################################################################################
#'
#' Prediction of the cluster memberships
#' 
#' This function gives the probabilities of classification for new observations by using the mixture model fit with the function \code{\link{VarSelCluster}}.
#' 
#' @param object instance of \code{\linkS4class{VSLCMresults}}.
#' @param newdata data.frame of the observations to classify.
#' @param type the type of prediction: probability of classification (probability) or the partition (partition)
#' 
#' @return Returns a matrix of the probabilities of classification.
#' 
#' @name predict
#' @rdname predict-methods
#' @docType methods
#' @exportMethod predict
#' @aliases predict predict,VSLCMresults-method
setMethod(f="predict",
          signature = c("VSLCMresults"),
          definition = function(object, newdata, type="probability"){
            #### Tests on the input arguments
            if (!(type %in% c("probability", "partition")))
              stop("type must be probability or partition")
            if ((ncol(newdata)!=object@data@d) || any(colnames(newdata) %in% object@data@var.names == FALSE) )
              stop("newdata must be contain the same features that the data used to fit the model")
            ####
            out <- proba.post(object, newdata)            
            colnames(out) <- paste("class", 1:ncol(out), sep = "-")
            if (type=="partition") out <- apply(out, 1, which.max)
            out
          } 
)