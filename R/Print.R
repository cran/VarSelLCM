########################################################################################################################
## Surcharge de la fonction print 
########################################################################################################################
#'
#' Print function.
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
    cat("Data set:\n   Number of individuals:", x@data@n,"\n")
    if (x@data@withContinuous){
      cat("   Number of continuous variables:", x@data@dataContinuous@d, "\n")
      val <- round(100*(1-mean(x@data@dataContinuous@notNA)),2)
      if (val>0)
        cat("   Percentile of missing values for the continuous variables:", val,"\n")
    }
    if (x@data@withInteger){
      cat("   Number of count variables:", x@data@dataInteger@d, "\n")
      val <- round(100*mean(1-x@data@dataInteger@notNA),2)
      if (val>0)
        cat("   Percentile of missing values for the integer variables:", val,"\n")
    }
    
    if (x@data@withCategorical){
      cat("   Number of categorical variables:", x@data@dataCategorical@d, "\n")
      miss <- 100*sum(sweep(is.na(x@data@dataCategorical@data),1,x@data@dataCategorical@weightdata,"*")) / (x@data@n * x@data@dataCategorical@d)
      if (miss>0)
        cat("   Percentile of missing values for the categorical variables:", miss,"\n")
    }
    cat("\n")
    
    cat("Model:\n   Number of components:", x@model@g, "\n")
    cat("   Model selection has been performed according to the", x@strategy@crit.varsel, " criterion \n")
    
    if (x@strategy@vbleSelec) cat("   Variable selection has been performed,", sum(x@model@omega)," (", round(100*sum(x@model@omega)/length(x@model@omega), 2),"% ) of the variables are relevant for clustering \n   \n")
    
    
    
    
    if ((length(x@criteria@degeneracyrate)==1)&&(x@criteria@degeneracyrate != 1)){
      cat("Information Criteria:\n")
      cat("   loglike:", x@criteria@loglikelihood,"\n")    
      cat("   AIC:    ", x@criteria@AIC,"\n")     
      cat("   BIC:    ", x@criteria@BIC,"\n")    
      cat("   ICL:    ", x@criteria@ICL,"\n")
      if ((x@strategy@crit.varsel=="MICL")&&(x@strategy@vbleSelec==TRUE)){        
        cat("   MICL:   ", x@criteria@MICL,"\n")    
        cat("   Best values has been found ", x@criteria@cvrate, "times\n")
      }  
    }
    cat("\n")
    if ((length(x@criteria@degeneracyrate)==1) && (x@criteria@degeneracyrate>0.1))
      cat("Warnings:\n  The rate of degeneracy for the EM algorithm is", x@criteria@degeneracyrate,"\n" )
  }
)