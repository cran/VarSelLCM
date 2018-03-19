########################################################################################################################
## Surcharge de la fonction summary 
########################################################################################################################
#'
#' Summary function.
#' 
#' This function gives the summary of an instance of  \code{\linkS4class{VSLCMresults}}.
#' 
#' @param object instance of \code{\linkS4class{VSLCMresults}}.
#' 
#' @name summary
#' @rdname summary-methods
#' @docType methods
#' @exportMethod summary
#' @aliases summary summary,VSLCMresults-method
setMethod(
  f="summary",
  signature = c("VSLCMresults"),
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
    
    cat("Model:\n   Number of components:", object@model@g, "\n")
    cat("   Model selection has been performed according to the", object@strategy@crit.varsel, " criterion \n")
    
    if (object@strategy@vbleSelec) cat("   Variable selection has been performed,", sum(object@model@omega)," (", round(100*sum(object@model@omega)/length(object@model@omega), 2),"% ) of the variables are relevant for clustering \n   \n")
    
     
    

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