ImputCont <- function(data, tik, param, method){
  output <- as.data.frame(data@data)
  for (j in 1:data@d){
    if (sum(data@notNA[,j])<data@n){
      who <- which(data@notNA[,j]==0)
      output[who, j] <- 0
      if (method=="postmean"){
        for (k in 1:ncol(tik)) output[who, j] <- output[who, j] + tik[who,k]*param@mu[j,k]
      }else if (method=="sampling"){
        z <- sample(1:ncol(tik), 1, prob=tik[who,])
        output[who, j] <- rnorm(1, param@mu[j,z], param@sd[j,z])
      }else{
        stop("argument method must be equal to postmean or sampling")
      }
    }  
  }
  return(output)
}

ImputInte <- function(data, tik, param, method){
  output <- as.data.frame(data@data)
  for (j in 1:data@d){
    if (sum(data@notNA[,j])<data@n){
      who <- which(data@notNA[,j]==0)
      output[who, j] <- 0
      if (method=="postmean"){
        for (k in 1:ncol(tik))  output[who, j] <- output[who, j] + tik[who,k]*param@lambda[j,k]
      }else if (method=="sampling"){
        z <- sample(1:ncol(tik), 1, prob=tik[who,])
        output[who, j] <- rpois(1, param@lambda[j,z])
      }else{
        stop("argument method must be equal to postmean or sampling")
      }
    }  
  }
  return(output)
}


ImputCate <- function(data, tik, param, method){
  output <- as.data.frame(data@data)
  for (j in 1:data@d){
    if (any(is.na(data@data[,j]))){
      who <- which(is.na(data@data[,j])==T)
      if (method=="postmean"){
        output[who, j] <- data@modalitynames[[j]][apply(tik[who,]%*%param@alpha[[j]],1,which.max)]
      }else if (method=="sampling"){
        z <- sample(1:ncol(tik), 1, prob=tik[who,])
        output[who, j] <- sample(data@modalitynames[[j]], 1, prob=param@alpha[[j]][z,])
      }else{
        stop("argument method must be equal to postmean or sampling")
      }
    }
  }
  return(output)
}

###################################################################################
###################################################################################
##' Imputation of missing values
##'
##' @description  
##' This function permits imputation of missing values in a dataset by using mixture model.
##' Two methods can be used for imputation:
##' \itemize{
#'  \item{posterior mean (method="postmean")}
#'  \item{sampling from the full conditionnal distribution (method="sampling")}
#'  }
##' 
##' @param obj an instance of \linkS4class{VSLCMresults}  which defines the model used for imputation.
##' @param newdata data.frame Dataset containing the missing values to impute.
##' @param method character definiting the method of imputation: "postmean" or "sampling"
##' 
##' @examples
##' # Data loading
##' data("heart")
##' 
##' # Clustering en 2 classes
##' results <- VarSelCluster(heart[,-13], 2)
##' 
##' # Data where missing values will be imputed
##' newdata <- heart[1:2,-13]
##' newdata[1,1] <- NA
##' newdata[2,2] <- NA
##' 
##' # Imputation
##' VarSelImputation(results, newdata)
##' 
##' @export
##'
##'
VarSelImputation <- function(obj, newdata, method="postmean"){
  #### Tests on the input arguments
  check.results(obj)        
  if (!(method %in% c("postmean", "sampling")))
    stop("method must be postmean or sampling")
  ####
  tik <- predict(obj, newdata)
  
  if (method=="postmean"){
    for (nom in colnames(newdata)){
      loc <- which(colnames(newdata)==nom)
      if (any(is.na(newdata[,loc]))){
        where <- which(is.na(newdata[,loc]))
        if (nom %in% rownames(obj@param@paramContinuous@mu)){
          who <- which(nom == rownames(obj@param@paramContinuous@mu))
          newdata[where, loc] <- as.numeric(tik[where, , drop=FALSE] %*% as.numeric(obj@param@paramContinuous@mu[who, , drop=FALSE]))
        }else if (nom %in% rownames(obj@param@paramInteger@lambda)){
          who <- which(nom == rownames(obj@param@paramInteger@lambda))
          newdata[where, loc] <- as.numeric(tik[where, , drop=FALSE] %*% as.numeric(obj@param@paramInteger@lambda[who, , drop=FALSE]))
        }else if (nom %in% names(obj@param@paramCategorical@alpha)){
          who <- which(nom ==  names(obj@param@paramCategorical@alpha))
          newdata[where, loc] <- obj@data@dataCategorical@modalitynames[[who]][apply(tik[where,]%*%obj@param@paramCategorical@alpha[[who]],1,which.max)]
        }
      }
    }
  }else{
    for (i in which(rowSums(is.na(newdata))>0)){
      zi <- sample(1:ncol(tik), 1, prob = tik[i,])
      for (nom in colnames(newdata)){
        loc <- which(colnames(newdata)==nom)
        if (is.na(newdata[i,loc])){
          if (nom %in% rownames(obj@param@paramContinuous@mu)){
            who <- which(nom == rownames(obj@param@paramContinuous@mu))
            newdata[i, loc] <- rnorm(1, obj@param@paramContinuous@mu[who, zi], obj@param@paramContinuous@sd[who, zi])
          }else if (nom %in% rownames(obj@param@paramInteger@lambda)){
            who <- which(nom == rownames(obj@param@paramInteger@lambda))
            newdata[i, loc] <- rpois(1, obj@param@paramInteger@lambda[who,zi])
          }else if (nom %in% names(obj@param@paramCategorical@alpha)){
            who <- which(nom ==  names(obj@param@paramCategorical@alpha))
            newdata[i, loc] <- sample(obj@data@dataCategorical@modalitynames[[who]], 1, prob = obj@param@paramCategorical@alpha[[who]][zi,])
          }
        }
      }
    }
  }
  return(newdata)         
}

