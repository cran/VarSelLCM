#' Constructor of \code{\linkS4class{VSLCMparamContinuous}} class
#' 
#' \describe{
#'   \item{pi}{numeric. Proportions of the mixture components.}
#'   \item{mu}{matrix. Mean for each component (column) and each variable (row).}
#'   \item{sd}{matrix. Standard deviation for each component (column) and each variable (row).}
#' }
#' @examples
#'   getSlots("VSLCMparamContinuous")
#' 
#' @name VSLCMparamContinuous-class
#' @rdname VSLCMparamContinuous-class
#' @exportClass VSLCMparamContinuous
setClass(
  Class = "VSLCMparamContinuous", 
  representation = representation(pi="numeric", mu="matrix", sd="matrix"), 
  prototype = prototype(pi=numeric(), mu=matrix(), sd=matrix())
)
##' Constructor of \code{\linkS4class{VSLCMparamInteger}} class
##' 
#' 
#' \describe{
#'   \item{pi}{numeric. Proportions of the mixture components.}
#'   \item{lambda}{matrix. Mean for each component (column) and each variable (row).}
#' }
#' @examples
#'   getSlots("VSLCMparamInteger")
#' 
#' @name VSLCMparamInteger-class
#' @rdname VSLCMparamInteger-class
#' @exportClass VSLCMparamInteger
setClass(
  Class = "VSLCMparamInteger", 
  representation = representation(pi="numeric", lambda="matrix"), 
  prototype = prototype(pi=numeric(), lambda=matrix())
)
##' Constructor of \code{\linkS4class{VSLCMparamCategorical}} class
##' 
#' \describe{
#'   \item{pi}{numeric. Proportions of the mixture components.}
#'   \item{alpha}{list. Parameters of the multinomial distributions.}
#' }
#' @examples
#'   getSlots("VSLCMparamCategorical")
#' 
#' @name VSLCMparamCategorical-class
#' @rdname VSLCMparamCategorical-class
#' @exportClass VSLCMparamCategorical
setClass(
  Class = "VSLCMparamCategorical", 
  representation = representation(pi="numeric", alpha="list"), 
  prototype = prototype(pi=numeric(), alpha=list())
)
########################################################################################################################
## Classe S4 VSLCMparamMixed contenant les parametres continus et categoriels
########################################################################################################################
###################################################################################
##' Constructor of \code{\linkS4class{VSLCMparam}} class
##'
##'  
##' \describe{
##'   \item{pi}{numeric. Proportions of the mixture components.}
##'   \item{paramContinuous}{\linkS4class{VSLCMparamContinuous}. Parameters of the continuous variables.}
##'   \item{paramInteger}{\linkS4class{VSLCMparamInteger}. Parameters of the integer variables.}
##'   \item{paramCategorical}{\linkS4class{VSLCMparamCategorical}. Parameters of the categorical variables.}
##' }
##'
##' @examples
##'   getSlots("VSLCMparam")
##'
##' @name VSLCMparam-class
##' @rdname VSLCMparam-class
##' @exportClass VSLCMparam
setClass(
  Class = "VSLCMparam", 
  representation = representation(pi="numeric", paramContinuous="VSLCMparamContinuous", paramInteger="VSLCMparamInteger", paramCategorical="VSLCMparamCategorical"), 
  prototype = prototype(pi=numeric(), paramContinuous=new("VSLCMparamContinuous"), paramInteger=new("VSLCMparamInteger"), paramCategorical=new("VSLCMparamCategorical"))
)
setGeneric(name= "convertparam", def = function(x) standardGeneric("convertparam"))

setMethod( f = "convertparam",
           signature(x="VSLCMparamCategorical"),
           definition = function(x) new("VSLCMparam",
                                        pi=x@pi,
                                        paramContinuous=new("VSLCMparamContinuous"), paramInteger=new("VSLCMparamInteger"), paramCategorical=x)
)


setMethod( f = "convertparam",
           signature(x="VSLCMparamContinuous"),
           definition = function(x) new("VSLCMparam",
                                        pi=x@pi,
                                        paramContinuous=x, paramInteger=new("VSLCMparamInteger"), paramCategorical=new("VSLCMparamCategorical"))
)


setMethod( f = "convertparam",
           signature(x="VSLCMparamInteger"),
           definition = function(x) new("VSLCMparam",
                                        pi=x@pi,
                                        paramContinuous=new("VSLCMparamContinuous"), paramInteger=x, paramCategorical=new("VSLCMparamCategorical"))
)

setMethod( f = "convertparam", signature(x="VSLCMparam"),  definition = function(x) x)

