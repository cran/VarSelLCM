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
    cat("Model:\n   Number of components:", object@model@g, "\n")
    cat("   Model selection has been performed according to the", object@strategy@crit.varsel, " criterion \n")
    
    if (object@strategy@vbleSelec) cat("   Variable selection has been performed,", sum(object@model@omega)," (", round(100*sum(object@model@omega)/length(object@model@omega), 2),"% ) of the variables are relevant for clustering \n   \n")
    
     

    cat("\n")
    if ((length(object@criteria@degeneracyrate)==1) && (object@criteria@degeneracyrate>0.1))
      cat("Warnings:\n  The rate of degeneracy for the EM algorithm is", object@criteria@degeneracyrate,"\n" )
  }
)