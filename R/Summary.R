setMethod(
  f="summary",
  signature = c("VSLCMresults"),
  definition = function(object){
    
    cat("Data set:\n   Number of individuals:", nrow(object@data), "\n   Number of variables:", ncol(object@data), "\n\n")
    
    cat("Model:\n   Number of components", object@model@g, "\n   Number of relevant variables for the clustering",sum(object@model@omega),"\n")
    if (sum(object@model@omega)>0)
      cat("\nNames of the relevant variables for the clustering:\n  ", colnames(object@data)[which(object@model@omega==1)])
    
    cat("\n\nInformation Criteria:\n")
    cat("   loglike:", object@criteria@likelihood,"\n")    
    cat("   BIC:    ", object@criteria@BIC,"\n")    
    cat("   ICL:    ", object@criteria@ICL,"\n")    
    cat("   MICL:   ", object@criteria@MICL,"\n")       
    
  }
)

