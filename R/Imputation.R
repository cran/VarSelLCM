ImputCont <- function(data, tik, param){
  output <- as.data.frame(data@data)
  for (j in 1:data@d){
    if (sum(data@notNA[,j])<data@n){
      who <- which(data@notNA[,j]==0)
      output[who, j] <- 0
      for (k in 1:ncol(tik))
        output[who, j] <- output[who, j] + tik[who,k]*param@mu[j,k]
    }  
  }
  return(output)
}

ImputInte <- function(data, tik, param){
  output <- as.data.frame(data@data)
  for (j in 1:data@d){
    if (sum(data@notNA[,j])<data@n){
      who <- which(data@notNA[,j]==0)
      output[who, j] <- 0
      for (k in 1:ncol(tik))
        output[who, j] <- output[who, j] + tik[who,k]*param@lambda[j,k]
    }  
  }
  return(output)
}


ImputCate <- function(data, tik, param){
  output <- as.data.frame(data@data)
  for (j in 1:data@d){
    if (any(is.na(data@data[,j]))){
      who <- which(is.na(data@data[,j])==T)
      output[who, j] <- data@modalitynames[[j]][apply(tik[who,]%*%param@alpha[[j]],1,which.max)]
    }
   # output[,j] <- as.factor(data@modalitynames[[j]][output[,j]])
  }
  return(output)
}

VarSelImputation <- function(obj, ind){
  if (missing(ind))
    ind <- 1:obj@data@n
  
  if (any( (ind %in% 1:obj@data@n)==FALSE))
    stop("Indice of individual is not correct!")
  
  output <- NULL
  
  if (class(obj)=="VSLCMresultsContinuous")
    output <- ImputCont(obj@data, obj@partitions@tik, obj@param)[ind,]
  else if (class(obj)=="VSLCMresultsInteger")
    output <- ImputInte(obj@data, obj@partitions@tik, obj@param)[ind,]
  else if (class(obj)=="VSLCMresultsCategorical")
    output <- ImputCate(obj@data, obj@partitions@tik, obj@param)[ind,]
  else if (class(obj)=="VSLCMresultsMixed"){
    output <- matrix(NA, length(ind), obj@data@d)
    colnames(output) <- names(obj@model@omega)
    output <- as.data.frame(output)
    col <- 1
    if (obj@data@withContinuous){
      output[,c(col:(col+obj@data@dataContinuous@d-1))] <- ImputCont(obj@data@dataContinuous, obj@partitions@tik, obj@param@paramContinuous)[ind,]
      col <- col + obj@data@dataContinuous@d
    }
    if (obj@data@withInteger){
      output[,c(col: (col+obj@data@dataInteger@d-1))] <- ImputInte(obj@data@dataInteger, obj@partitions@tik, obj@param@paramInteger)[ind,]
      col <- col + obj@data@dataInteger@d
    }
    if (obj@data@withCategorical){
      output[,c(col: (col+obj@data@dataCategorical@d-1))] <- ImputCate(obj@data@dataCategorical, obj@partitions@tik, obj@param@paramCategorical)[ind,]
      col <- col + obj@data@dataCategorical@d
    }
  }  else
    stop("obj don't arise from function VarSelCluster")
  
  or <- 1:obj@data@d
  for (j in 1:obj@data@d){
    or[j] <- which(colnames(output)==obj@data@var.names[j])
  }
  output <- output[,or]
  
  return(output)         
}

