setMethod(
  f="print",
  signature = c("VSLCMresults"),
  definition = function(x){
    
    summary(x)
    
    cat("\n Parameters per class:")
    for (k in 1:x@model@g){
      if (k>1){
        cat("*******************\n")
      }
      cat("Class",k,"\n")
      cat("Proportion:",x@parameters@proportions[k],"\n")
      tmp <- data.frame(mean=x@parameters@means[k,], variance=x@parameters@variances[k,])
      print(tmp)
      cat("\n")
    }
    
  }
)

