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
##' Imputation function based on the mixture model. Two methods can be used: missing values can be 
##' imputed by their posterior mean (args "postmean") or by a sampling from their full conditionnal
##' distribution (args "sampling").
##' 
##' @param obj an instance of \linkS4class{VSLCMresults} returned by function \link{VarSelCluster}.
##' @param method character definiting the method of imputation: "postmean" or "sampling"
##' 
##' @examples
##' \dontrun{
##' # Data loading
##' data("heart")
##' # Clustering en 2 classes
##' heart[1,1] <- NA
##' results <- VarSelCluster(heart[,-13], 2)
##' # Opening Shiny application to easily see the results
##' VarSelImputation(results)[1,1]
##' }
##' 
##' @export
##'
##'
VarSelImputation <- function(obj, method="postmean"){
  output <- NULL
  if (class(obj)=="VSLCMresults"){
    output <- matrix(NA, obj@data@n, obj@data@d)
    colnames(output) <- names(obj@model@omega)
    output <- as.data.frame(output)
    col <- 1
    if (obj@data@withContinuous){
      output[,c(col:(col+obj@data@dataContinuous@d-1))] <- ImputCont(obj@data@dataContinuous, obj@partitions@tik, obj@param@paramContinuous, method)
      col <- col + obj@data@dataContinuous@d
    }
    if (obj@data@withInteger){
      output[,c(col: (col+obj@data@dataInteger@d-1))] <- ImputInte(obj@data@dataInteger, obj@partitions@tik, obj@param@paramInteger, method)
      col <- col + obj@data@dataInteger@d
    }
    if (obj@data@withCategorical){
      output[,c(col: (col+obj@data@dataCategorical@d-1))] <- ImputCate(obj@data@dataCategorical, obj@partitions@tik, obj@param@paramCategorical, method)
      col <- col + obj@data@dataCategorical@d
    }
  } else
    stop("obj doesn't arise from function VarSelCluster")
  
  or <- 1:obj@data@d
  for (j in 1:obj@data@d){
    or[j] <- which(colnames(output)==obj@data@var.names[j])
  }
  output <- output[,or]
  
  return(output)         
}

