########################################################################################################################
## Surcharge de la fonction summary pour les objects de classe S4 VSLCMresultsContinuous et VSLCMresultsCategorical
########################################################################################################################
#'
#' Summary function.
#' 
#' This function gives the summary of an instance of \code{\linkS4class{VSLCMresultsContinuous}}, \code{\linkS4class{VSLCMresultsInteger}}, \code{\linkS4class{VSLCMresultsCategorical}} or \code{\linkS4class{VSLCMresultsMixed}}.
#' 
#' @param object instance of  \code{\linkS4class{VSLCMresultsContinuous}}, \code{\linkS4class{VSLCMresultsInteger}}, \code{\linkS4class{VSLCMresultsCategorical}} or \code{\linkS4class{VSLCMresultsMixed}}..
#' 
#' @name summary
#' @rdname summary-methods
#' @docType methods
#' @exportMethod summary
## Surcharge pour VSLCMresultsContinuous
#' @rdname summary-methods
#' @aliases summary summary,VSLCMresultsContinuous-method
setMethod(
  f="summary",
  signature = c("VSLCMresultsContinuous"),
  definition = function(object){
    cat("Data set:\n   Number of individuals:", object@data@n,"\n")
    cat("   Number of continuous variables:", object@data@d, "\n")
    val <- round(100*(1-mean(object@data@notNA)),2)
    if (val>0) cat("   Percentile of missing values:", val,"\n\n")
    cat("Model:\n   Number of components:", object@model@g, "\n   Number of relevant variables for the clustering:",sum(object@model@omega)," (", 100*sum(object@model@omega)/length(object@model@omega),"% )")
    if (object@strategy@vbleSelec){
      cat("\n   The variable selection has been performed according to the", object@strategy@crit.varsel, " criterion \n")
    }else{
      cat("\n   No variable selection has been performed \n")
    }
    if (sum(object@model@omega)>0){
      if (length(length(object@model@names.relevant))>6){
        cat("   Names of the first six relevant variables for the clustering: ")        
      }else{
        cat("   Names of the relevant variables for the clustering: ") 
      }
      cat(object@model@names.relevant[1:min(6,length(object@model@names.relevant))], "\n\n")
    }   
    if ((length(object@criteria@degeneracyrate)==1)&&(object@criteria@degeneracyrate != 1)){
      cat("Information Criteria:\n")
      cat("   loglike:", object@criteria@loglikelihood,"\n")    
      cat("   AIC:    ", object@criteria@AIC,"\n")    
      cat("   BIC:    ", object@criteria@BIC,"\n")    
      cat("   ICL:    ", object@criteria@ICL,"\n") 
      if (object@strategy@vbleSelec){
        
        cat("   MICL:   ", object@criteria@MICL,"\n")    
        cat("   Best values has been found ", object@criteria@cvrate, "times\n")
      }  
    }
    cat("\n")
    if ((length(object@criteria@degeneracyrate)==1) && (object@criteria@degeneracyrate>0.1))
      cat("Warnings:\n  The rate of degeneracy for the EM algorithm is", object@criteria@degeneracyrate,"\n" )
    
  }
)


## Surcharge pour VSLCMresultsInteger
#' @rdname summary-methods
#' @aliases summary summary,VSLCMresultsInteger-method
setMethod(
  f="summary",
  signature = c("VSLCMresultsInteger"),
  definition = function(object){
    cat("Data set:\n   Number of individuals:", object@data@n,"\n")
    cat("   Number of count variables:", object@data@d, "\n")
    val <- round(100*(1-mean(object@data@notNA)),2)
    if (val>0)
      cat("   Percentile of missing values:", ,"\n\n")
    cat("Model:\n   Number of components:", object@model@g, "\n   Number of relevant variables for the clustering:",sum(object@model@omega)," (", 100*sum(object@model@omega)/length(object@model@omega),"% ) \n")
    if (object@strategy@vbleSelec){
      cat("\n   The variable selection has been performed according to the", object@strategy@crit.varsel, " criterion \n")
    }else{
      cat("\n   No variable selection has been performed \n")
    }
    if (sum(object@model@omega)>0){
      if (length(length(object@model@names.relevant))>6){
        cat("   Names of the first six relevant variables for the clustering: ")        
      }else{
        cat("   Names of the relevant variables for the clustering: ") 
      }
      cat(object@model@names.relevant[1:min(6,length(object@model@names.relevant))], "\n\n")
    }     
    if ((length(object@criteria@degeneracyrate)==1)&&(object@criteria@degeneracyrate != 1)){
      cat("Information Criteria:\n")
      cat("   loglike:", object@criteria@loglikelihood,"\n")    
      cat("   AIC:    ", object@criteria@AIC,"\n")    
      cat("   BIC:    ", object@criteria@BIC,"\n")    
      cat("   ICL:    ", object@criteria@ICL,"\n") 
      if (object@strategy@vbleSelec){
        
        cat("   MICL:   ", object@criteria@MICL,"\n")    
        cat("   Best values has been found ", object@criteria@cvrate, "times\n")
      } 
    }
    cat("\n")
    if ((length(object@criteria@degeneracyrate)==1) && (object@criteria@degeneracyrate>0.1))
      cat("Warnings:\n  The rate of degeneracy for the EM algorithm is", object@criteria@degeneracyrate,"\n" )
    
  }
)

## Surcharge pour VSLCMresultsCategorical
#' @rdname summary-methods
#' @aliases summary summary,VSLCMresultsCategorical-method
setMethod(
  f="summary",
  signature = c("VSLCMresultsCategorical"),
  definition = function(object){
    cat("Data set:\n   Number of individuals:", object@data@n,"\n")
    cat("   Number of profiles:", nrow(object@data@shortdata),"\n")
    cat("   Number of categorical variables:", object@data@d, "\n")
    miss <- sum(sweep(is.na(object@data@shortdata),1,object@data@weightdata,"*")) / (object@data@n * object@data@d)
    if (miss>0)
      cat("   Percentile of missing values:", miss,"\n\n")
    cat("Model:\n   Number of components:", object@model@g, "\n   Number of relevant variables for the clustering:",sum(object@model@omega)," (", 100*sum(object@model@omega)/length(object@model@omega),"% ) \n")
    if (object@strategy@vbleSelec){
      cat("\n   The variable selection has been performed according to the", object@strategy@crit.varsel, " criterion \n")
    }else{
      cat("\n   No variable selection has been performed \n")
    }
    if (sum(object@model@omega)>0){
      if (length(length(object@model@names.relevant))>6){
        cat("   Names of the first six relevant variables for the clustering: ")        
      }else{
        cat("   Names of the relevant variables for the clustering: ") 
      }
      cat(object@model@names.relevant[1:min(6,length(object@model@names.relevant))], "\n\n")
    }    
    cat("Information Criteria:\n")
    cat("   loglike:", object@criteria@loglikelihood,"\n")    
    cat("   AIC:    ", object@criteria@AIC,"\n")     
    cat("   BIC:    ", object@criteria@BIC,"\n")    
    cat("   ICL:    ", object@criteria@ICL,"\n")
    if (object@strategy@vbleSelec){
      
      cat("   MICL:   ", object@criteria@MICL,"\n")    
      cat("   Best values has been found ", object@criteria@cvrate, "times\n")
    }  
  }
)

## Surcharge pour VSLCMresultsMixed
#' @rdname summary-methods
#' @aliases summary summary,VSLCMresultsMixed-method
setMethod(
  f="summary",
  signature = c("VSLCMresultsMixed"),
  definition = function(object){
    cat("Data set:\n   Number of individuals:", object@data@n,"\n")
    if (object@data@withContinuous){
      cat("   Number of continuous variables:", object@data@dataContinuous@d, "\n")
      val <- round(100*(1-mean(object@data@dataContinuous@notNA)),2)
      if (val>0)
        cat("   Percentile of missing values for the continuous variables:", val,"\n")
    }
    if (object@data@withInteger){
      cat("   Number of count variables:", object@data@dataInteger@d, "\n")
      val <- round(100*mean(1-object@data@dataInteger@notNA),2)
      if (val>0)
        cat("   Percentile of missing values for the integer variables:", val,"\n")
    }
    
    if (object@data@withCategorical){
      cat("   Number of categorical variables:", object@data@dataCategorical@d, "\n")
      miss <- 100*sum(sweep(is.na(object@data@dataCategorical@data),1,object@data@dataCategorical@weightdata,"*")) / (object@data@n * object@data@dataCategorical@d)
      if (miss>0)
        cat("   Percentile of missing values for the categorical variables:", miss,"\n")
    }
    cat("\n")
    
    cat("Model:\n   Number of components:", object@model@g, "\n   Number of relevant variables for the clustering:",sum(object@model@omega)," (", 100*sum(object@model@omega)/length(object@model@omega),"% ) \n")
    if (object@strategy@vbleSelec){
      cat("   The variable selection has been performed according to the", object@strategy@crit.varsel, " criterion \n")
    }else{
      cat("   No variable selection has been performed \n")
    }
    if (sum(object@model@omega)>0){
      if (length(length(object@model@names.relevant))>6){
        cat("   Names of the first six relevant variables for the clustering: ")        
      }else{
        cat("   Names of the relevant variables for the clustering: ") 
      }
      cat(object@model@names.relevant[1:min(6,length(object@model@names.relevant))], "\n\n")
    }     
    if ((length(object@criteria@degeneracyrate)==1)&&(object@criteria@degeneracyrate != 1)){
      cat("Information Criteria:\n")
      cat("   loglike:", object@criteria@loglikelihood,"\n")    
      cat("   AIC:    ", object@criteria@AIC,"\n")     
      cat("   BIC:    ", object@criteria@BIC,"\n")    
      cat("   ICL:    ", object@criteria@ICL,"\n")
      if ((object@strategy@crit.varsel=="MICL")&&(object@strategy@vbleSelec==TRUE)){        
        cat("   MICL:   ", object@criteria@MICL,"\n")    
        cat("   Best values has been found ", object@criteria@cvrate, "times\n")
      }  
    }
    cat("\n")
    if ((length(object@criteria@degeneracyrate)==1) && (object@criteria@degeneracyrate>0.1))
      cat("Warnings:\n  The rate of degeneracy for the EM algorithm is", object@criteria@degeneracyrate,"\n" )
  }
)