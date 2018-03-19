###################################################################################
###################################################################################
##' Convert outputs of old versions of VarSelLCM
##'
##' @description  
##' This function is used to convert output obtained by old versions (before 2.1.0) if VarSelLCM
##' 
##' @param obj an instance returned by function \link{VarSelCluster}.
##' 
##' @examples
##' \dontrun{
##' ##' # Data loading
##' data("heart")
##' # Clustering en 2 classes
##' heart[1,1] <- NA
##' results <- VarSelCluster(heart[,-13], 2)
##' # Here, it is necessary only if results was provided by an older version of the VarSelLCM
##' update <- VarSelConvert(results)
##' }
##' 
##' @export
##'
##'
VarSelConvert <- function(obj){
  if (!(class(obj) %in% c("VSLCMresultsContinuous", "VSLCMresultsInteger", "VSLCMresultsCategorical", "VSLCMresultsMixed")))
    stop("This function is used to convert objects return by old version of the R package VarSelLCM, your input argument does not arise from an old version of this package")
  new("VSLCMresults", data=convertdata(obj@data), criteria=obj@criteria, partitions=obj@partitions, model=obj@model, strategy=obj@strategy, param=convertparam(obj@param))
}