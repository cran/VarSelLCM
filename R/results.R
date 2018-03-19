#' Constructor of \code{\linkS4class{VSLCMresultsContinuous}} class
#' 
#'  \describe{
#'   \item{data}{\linkS4class{VSLCMdataContinuous}. Results relied to the data.}
#'   \item{criteria}{\linkS4class{VSLCMcriteria}. Results relied to the information criteria.}
#'   \item{partitions}{\linkS4class{VSLCMpartitions}. Results relied to the partitions.}
#'   \item{model}{\linkS4class{VSLCMmodel}. Results relied to the selected model.}
#'   \item{strategy}{\linkS4class{VSLCMstrategy}. Results relied to the tune parameters.}
#'   \item{param}{\linkS4class{VSLCMparamContinuous}. Results relied to the parameters.}
#' }
#' @examples
#' getSlots("VSLCMresultsContinuous")
#' 
#' @name VSLCMresultsContinuous-class
#' @rdname VSLCMresultsContinuous-class
#' @exportClass VSLCMresultsContinuous
setClass(
  Class = "VSLCMresultsContinuous", 
  representation = representation(data="VSLCMdataContinuous", criteria="VSLCMcriteria", partitions="VSLCMpartitions",
                                  model="VSLCMmodel", strategy="VSLCMstrategy", param="VSLCMparamContinuous"), 
  prototype = prototype(data=new("VSLCMdataContinuous"), criteria=new("VSLCMcriteria"), partitions=new("VSLCMpartitions"),
                        model=new("VSLCMmodel"), strategy=new("VSLCMstrategy"), param=new("VSLCMparamContinuous"))
)

#' Constructor of \code{\linkS4class{VSLCMresultsInteger}} class
#' 
#' \describe{
#'   \item{data}{\linkS4class{VSLCMdataInteger}. Results relied to the data.}
#'   \item{criteria}{\linkS4class{VSLCMcriteria}. Results relied to the information criteria.}
#'   \item{partitions}{\linkS4class{VSLCMpartitions}. Results relied to the partitions.}
#'   \item{model}{\linkS4class{VSLCMmodel}. Results relied to the selected model.}
#'   \item{strategy}{\linkS4class{VSLCMstrategy}. Results relied to the tune parameters.}
#'   \item{param}{\linkS4class{VSLCMparamInteger}. Results relied to the parameters.}
#' }
#' @examples
#' getSlots("VSLCMresultsInteger")
#' 
#' @name VSLCMresultsInteger-class
#' @rdname VSLCMresultsInteger-class
#' @exportClass VSLCMresultsInteger
setClass(
  Class = "VSLCMresultsInteger", 
  representation = representation(data="VSLCMdataInteger", criteria="VSLCMcriteria", partitions="VSLCMpartitions",
                                  model="VSLCMmodel", strategy="VSLCMstrategy", param="VSLCMparamInteger"), 
  prototype = prototype(data=new("VSLCMdataInteger"), criteria=new("VSLCMcriteria"), partitions=new("VSLCMpartitions"),
                        model=new("VSLCMmodel"), strategy=new("VSLCMstrategy"), param=new("VSLCMparamInteger"))
)

#' Constructor of \code{\linkS4class{VSLCMresultsCategorical}} class
#' 
#'  \describe{
#'   \item{data}{\linkS4class{VSLCMdataCategorical}. Results relied to the data.}
#'   \item{criteria}{\linkS4class{VSLCMcriteria}. Results relied to the information criteria.}
#'   \item{partitions}{\linkS4class{VSLCMpartitions}. Results relied to the partitions.}
#'   \item{model}{\linkS4class{VSLCMmodel}. Results relied to the selected model.}
#'   \item{strategy}{\linkS4class{VSLCMstrategy}. Results relied to the tune parameters.}
#'   \item{param}{\linkS4class{VSLCMparamCategorical}. Results relied to the parameters.}
#' }
#' @examples
#'   getSlots("VSLCMresultsCategorical")
#' 
#' @name VSLCMresultsCategorical-class
#' @rdname VSLCMresultsCategorical-class
#' @exportClass VSLCMresultsCategorical
setClass(
  Class = "VSLCMresultsCategorical", 
  representation = representation(data="VSLCMdataCategorical", criteria="VSLCMcriteria", partitions="VSLCMpartitions",
                                  model="VSLCMmodel", strategy="VSLCMstrategy", param="VSLCMparamCategorical"), 
  prototype = prototype(data=new("VSLCMdataCategorical"), criteria=new("VSLCMcriteria"), partitions=new("VSLCMpartitions"),
                        model=new("VSLCMmodel"), strategy=new("VSLCMstrategy"), param=new("VSLCMparamCategorical"))
)




########################################################################################################################
## Classe S4 VSLCMresultsMixed
########################################################################################################################
###################################################################################
##' Constructor of \code{\linkS4class{VSLCMresults}} class
##'
##'  
##' \describe{
##'   \item{data}{\linkS4class{VSLCMdata}. Results relied to the data.}
##'   \item{criteria}{\linkS4class{VSLCMcriteria}. Results relied to the information criteria.}
##'   \item{partitions}{\linkS4class{VSLCMpartitions}. Results relied to the partitions.}
##'   \item{model}{\linkS4class{VSLCMmodel}. Results relied to the selected model.}
##'   \item{strategy}{\linkS4class{VSLCMstrategy}. Results relied to the tune parameters.}
##'   \item{param}{\linkS4class{VSLCMparam}. Results relied to the parameters.}
##' }
##'
#' @examples
#'   getSlots("VSLCMresults")
#'
#' @name VSLCMresults-class
#' @rdname VSLCMresults-class
#' @exportClass VSLCMresults
setClass(
  Class = "VSLCMresults", 
  representation = representation(data="VSLCMdata", criteria="VSLCMcriteria", partitions="VSLCMpartitions",
                                  model="VSLCMmodel", strategy="VSLCMstrategy", param="VSLCMparam"), 
  prototype = prototype(data=new("VSLCMdata"), criteria=new("VSLCMcriteria"), partitions=new("VSLCMpartitions"),
                        model=new("VSLCMmodel"), strategy=new("VSLCMstrategy"), param=new("VSLCMparam"))
)


