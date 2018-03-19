########################################################################################################################
## Surcharge de la fonction print 
########################################################################################################################
#'
#' Summary function.
#' 
#' This function gives the print of an instance of  \code{\linkS4class{VSLCMresults}}.
#' 
#' @param x instance of \code{\linkS4class{VSLCMresults}}.
#' 
#' @name print
#' @rdname print-methods
#' @docType methods
#' @exportMethod print
#' @aliases print print,VSLCMresults-method

setMethod(
  f="print",
  signature = c("VSLCMresults"),
  definition = function(x){
    summary(x)    
    if (x@criteria@degeneracyrate != 1){
      cat("**************************************\n")
      cat("\nParameters of the relevant variables:\n")
      for (k in 1:x@model@g){
        if (k>1){
          cat("*******************\n")
        }
        cat("Class",k,"\n")
        cat("Proportion:",x@param@pi[k],"\n")
        if (x@data@withContinuous){
          
          tmp <- data.frame(mean=x@param@paramContinuous@mu[,k], sd=x@param@paramContinuous@sd[,k])
          if (sum(x@model@omega[names(x@model@omega)%in%rownames(tmp)])>0){
            cat("Parameters of continuous variables \n")
            keep <- as.data.frame(tmp[which(x@model@omega[names(x@model@omega)%in%rownames(tmp)]==1),])
            rownames(keep) <- rownames(tmp)[which(x@model@omega[names(x@model@omega)%in%rownames(tmp)]==1)]
            colnames(keep) <- c("mean", "sd")
            print(keep)
          }
            
          cat("\n")
        }
        if (x@data@withInteger){
          tmp <- data.frame(lambda=x@param@paramInteger@lambda[,k])
          rownames(tmp)  <- rownames(x@param@paramInteger@lambda)
          if (sum(x@model@omega[names(x@model@omega)%in%rownames(tmp)])>0){
            cat("Parameters of count variables \n")
            keep <- as.data.frame(tmp[which(x@model@omega[names(x@model@omega)%in%rownames(tmp)]==1),])
            rownames(keep) <- rownames(tmp)[which(x@model@omega[names(x@model@omega)%in%rownames(tmp)]==1)]         
            colnames(keep) <- "lambda"
            print(keep)
          }
          cat("\n")
        }
        
        if (x@data@withCategorical){
          maxcol <- 0
          for (j in 1:x@data@dataCategorical@d) maxcol <- max(maxcol, length(x@data@dataCategorical@modalitynames[[j]]))
          alpha <- matrix(0,x@data@dataCategorical@d, maxcol)
          for (j in 1:x@data@dataCategorical@d)  alpha[j,1:length(x@data@dataCategorical@modalitynames[[j]])] <- round(x@param@paramCategorical@alpha[[j]][k,],6)
          for (j in 1:ncol(alpha)) alpha[,j] <- as.character(alpha[,j])
          for (j in 1:x@data@dataCategorical@d){
            if (length(x@data@dataCategorical@modalitynames[[j]])<maxcol)
              alpha[j, (length(x@data@dataCategorical@modalitynames[[j]])+1):maxcol] <- rep(".", maxcol-length(x@data@dataCategorical@modalitynames[[j]]))
            
          }
          alpha <- data.frame(alpha)
          colnames(alpha) <- paste("Level",1:ncol(alpha),sep=".")
          rownames(alpha) <- names(x@param@paramCategorical@alpha)
          if (sum(x@model@omega[names(x@model@omega)%in%rownames(alpha)])>0){
            cat("Parameters of categorical variables \n")      
            keep <- alpha[which(x@model@omega[names(x@model@omega)%in%rownames(alpha)]==1),]
            rownames(keep) <- rownames(alpha)[which(x@model@omega[names(x@model@omega)%in%rownames(alpha)]==1)]
            print(keep) 
          }
          cat("\n")
        }
        
        cat("\n")
      }
      if (any(x@model@omega==0)){
        cat("**************************************\n")
        cat("\nParameters of the irrelevant variables:\n")
        if (x@data@withContinuous){
          tmp <- data.frame(mean=x@param@paramContinuous@mu[,k], sd=x@param@paramContinuous@sd[,k])
          if (any(x@model@omega[names(x@model@omega)%in%rownames(tmp)]==0)){
            cat("Parameters of countinous variables \n")
            keep <- as.data.frame(tmp[which(x@model@omega[names(x@model@omega)%in%rownames(tmp)]==0),])
            rownames(keep) <- rownames(tmp)[which(x@model@omega[names(x@model@omega)%in%rownames(tmp)]==0)]
            colnames(keep) <- c("mean", "sd")
            print(keep)
          }
          
          cat("\n")
        }
        if (x@data@withInteger){
          tmp <- data.frame(lambda=x@param@paramInteger@lambda[,k])
          rownames(tmp)  <- rownames(x@param@paramInteger@lambda)
          if (any(x@model@omega[names(x@model@omega)%in%rownames(tmp)]==0)){
            cat("Parameters of count variables \n")
            keep <- as.data.frame(tmp[which(x@model@omega[names(x@model@omega)%in%rownames(tmp)]==0),])
            rownames(keep) <- rownames(tmp)[which(x@model@omega[names(x@model@omega)%in%rownames(tmp)]==0)]
            colnames(keep) <- "lambda"
            print(keep)
          }
          cat("\n")
        }
        
        if (x@data@withCategorical){ 
          maxcol <- 0
          for (j in 1:x@data@dataCategorical@d) maxcol <- max(maxcol, length(x@data@dataCategorical@modalitynames[[j]]))
          alpha <- matrix(0,x@data@dataCategorical@d, maxcol)
          for (j in 1:x@data@dataCategorical@d)  alpha[j,1:length(x@data@dataCategorical@modalitynames[[j]])] <- round(x@param@paramCategorical@alpha[[j]][k,],6)
          for (j in 1:ncol(alpha)) alpha[,j] <- as.character(alpha[,j])
          for (j in 1:x@data@dataCategorical@d){
            if (length(x@data@dataCategorical@modalitynames[[j]])<maxcol)
              alpha[j, (length(x@data@dataCategorical@modalitynames[[j]])+1):maxcol] <- rep(".", maxcol-length(x@data@dataCategorical@modalitynames[[j]]))
            
          }
          alpha <- data.frame(alpha)
          colnames(alpha) <- paste("Level",1:ncol(alpha),sep=".")
          rownames(alpha) <- names(x@param@paramCategorical@alpha)
          if (any(x@model@omega[names(x@model@omega)%in%rownames(alpha)]==0)){
            cat("Parameters of categorical variables \n")      
            keep <- alpha[which(x@model@omega[names(x@model@omega)%in%rownames(alpha)]==0),]
            rownames(keep) <- rownames(alpha)[which(x@model@omega[names(x@model@omega)%in%rownames(alpha)]==0)]
            print(keep) 
          }
          cat("\n")
        }
      }
      
      
    }
  }
)