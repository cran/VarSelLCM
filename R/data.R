
########################################################################################################################## S4 classes for the data################################################################################################################################################################################################################################################## S4 VSLCMdataContinuous: class for continuous data sets#######################################################################################################################Constructor of [\code{\linkS4class{VSLCMdataContinuous}}] class

#' Constructor of \code{\linkS4class{VSLCMdataContinuous}} class
#' 
#' \describe{ 
#'   \item{n}{number of observations}
#'   \item{d}{number of variables}
#'   \item{data}{matrix of observations (one row = one observation)}
#'   \item{notNA}{matrix of logical (1:observed, 0:unobserved)}
#'   \item{priors}{hyper-parameters of the prior distributions}
#' }
#' @examples
#'   getSlots("VSLCMdataContinuous")
#'   
#' 
#' @name VSLCMdataContinuous-class
#' @rdname VSLCMdataContinuous-class
#' @exportClass VSLCMdataContinuous



setClass(
  Class = "VSLCMdataContinuous", 
  representation = representation(
    n="numeric",
    d="numeric",
    data="matrix",
    notNA="matrix",
    priors="matrix"
  ), 
  prototype = prototype(
    n=numeric(),
    d=numeric(),
    data=matrix(),
    notNA=matrix(),
    priors=matrix()
  )
)
#' Constructor of \code{\linkS4class{VSLCMdataInteger}} class
#' 
#' \describe{
#'   \item{n}{number of observations}
#'   \item{d}{number of variables}
#'   \item{data}{matrix of observations (one row = one observation)}
#'   \item{notNA}{matrix of logical (1:observed, 0:unobserved)}
#'   \item{priors}{hyper-parameters of the prior distributions}
#' }
#'  @examples
#'   getSlots("VSLCMdataInteger")
#' 
#' @name VSLCMdataInteger-class
#' @rdname VSLCMdataInteger-class
#' @exportClass VSLCMdataInteger

setClass(
  Class = "VSLCMdataInteger", 
  representation = representation(
    n="numeric",
    d="numeric",
    data="matrix",
    notNA="matrix",
    priors="matrix"
  ), 
  prototype = prototype(
    n=numeric(),
    d=numeric(),
    data=matrix(),
    notNA=matrix(),
    priors=matrix()
  )
)
#' Constructor of \code{\linkS4class{VSLCMdataCategorical}} class
#' 
#' \describe{
#'   \item{n}{number of observations}
#'   \item{d}{number of variables}
#'   \item{data}{matrix of observations (one row = one observation)}
#'   \item{shortdata}{matrix of unique profils}
#'   \item{weightdata}{weights of profils}
#'   \item{modalitynames}{names of levels}
#' }
#' 
#' @examples
#'   getSlots("VSLCMdataCategorical")
#' 
#' @name VSLCMdataCategorical-class
#' @rdname VSLCMdataCategorical-class
#' @exportClass VSLCMdataCategorical
setClass(
  Class = "VSLCMdataCategorical", 
  representation = representation(
    n="numeric",
    d="numeric",
    data="matrix",
    shortdata="matrix",
    weightdata="numeric",
    modalitynames="list"
  ), 
  prototype = prototype(
    n=numeric(),
    d=numeric(),
    data=matrix(),
    shortdata=matrix(),
    weightdata=numeric(),
    modalitynames=list()
  )
)



########################################################################################################################
########################################################################################################################
## VSLCdata: class of  data set
########################################################################################################################
##' Constructor of \code{\linkS4class{VSLCMdata}} class
##' 
##' \describe{
##'   \item{n}{number of observations}
##'   \item{d}{number of variables}
##' \item{withContinuous}{logical indicating if some variables are continuous}
##' \item{withInteger}{logical indicating if some variables are integer}
##' \item{withCategorica}{logical indicating if some variables are categorical} 
##' \item{dataContinuous}{instance of VSLCMdataContinuous containing the continuous data}
##' \item{dataInteger}{instance of VSLCMdataContinuous containing the integer data}
##' \item{dataCategorical}{instance of VSLCMdataContinuous containing the categorical data}
##' \item{var.names}{labels of the variables} 
##' }
##'
#' @examples
#'   getSlots("VSLCMdata")
#'
#' @name VSLCMdata-class
#' @rdname VSLCMdata-class
#' @exportClass VSLCMdata
setClass(
  Class = "VSLCMdata", 
  representation = representation(
    n="numeric",
    d="numeric",
    withContinuous="logical",
    withInteger="logical",
    withCategorical="logical",
    dataContinuous="VSLCMdataContinuous",
    dataInteger="VSLCMdataInteger",
    dataCategorical="VSLCMdataCategorical",
    var.names="character"
  ), 
  prototype = prototype(
    n=numeric(),
    d=numeric(),
    withContinuous=logical(),
    withInteger=logical(),
    withCategorical=logical(),
    dataContinuous=new("VSLCMdataContinuous"),
    dataInteger=new("VSLCMdataInteger"),
    dataCategorical=new("VSLCMdataCategorical"),
    var.names=character()
  )
)



setGeneric(name= "convertdata", def = function(x) standardGeneric("convertdata"))

setMethod( f = "convertdata",
           signature(x="VSLCMdataCategorical"),
           definition = function(x) new("VSLCMdata", n=x@n, d=x@d, 
                                        withContinuous=FALSE,  withInteger=FALSE, withCategorical=TRUE,
                                        dataContinuous=new("VSLCMdataContinuous"), dataInteger=new("VSLCMdataInteger"), dataCategorical=x, var.names=colnames(x@data))
)

setMethod( f = "convertdata",
           signature(x="VSLCMdataContinuous"),
           definition = function(x) new("VSLCMdata", n=x@n, d=x@d, 
                                        withContinuous=TRUE,  withInteger=FALSE, withCategorical=FALSE,
                                        dataContinuous=x, dataInteger=new("VSLCMdataInteger"), dataCategorical=new("VSLCMdataCategorical"), var.names=colnames(x@data))
)

setMethod( f = "convertdata",
           signature(x="VSLCMdataInteger"),
           definition = function(x) new("VSLCMdata", n=x@n, d=x@d, 
                                        withContinuous=FALSE,  withInteger=TRUE, withCategorical=FALSE,
                                        dataContinuous=new("VSLCMdataContinuous"), dataInteger=x, dataCategorical=new("VSLCMdataCategorical"), var.names=colnames(x@data))
)

setMethod( f = "convertdata",
           signature(x="VSLCMdata"),
           definition = function(x) x)

########################################################################################################################
## La fonction VSLCMdata permet de construire un objet de class S4 VSLCMdataContinuous ou VSLCMdataCategorical en fonction
## de la nature des variables
########################################################################################################################
VSLCMdata <- function(x, redquali=TRUE){
  # Ajout d'un nom de variable si celui-ci est manquant
  if (is.null(colnames(x))) colnames(x) <- paste("X",1:ncol(x), sep="")
  if (max(table(colnames(x)))>1) stop("At least two variables have the same name!")
  n <- nrow(x)
  d <- ncol(x)
  # recherche des indices de variables numeric et factor
  type <- numeric()
  for (j in 1:d) type[j] <- class(x[,j])
  idxcont <- which(type=="numeric")
  idxinte <- which(type=="integer")
  idxcat <- which(type=="factor")
  if ((all(type %in% c("numeric", "integer", "factor"))==FALSE))
    stop("At least one variable is neither numeric, integer nor factor!")
  
  # cas des variables categorielles
  if ( (length(idxcat) == d) ){
    shortdata <- matrix(NA, n, d)
    for (j in 1:d){
      lev <- levels(x[,j])
      repere <- 0
      for (h in 1:length(lev)){
        repere <- repere + 1
        shortdata[which(x[,j]==lev[h]),j] <- repere
      }
    }
    weightdata <- rep(1, n)
    ## Pour travailler avec Armadillo on rempli artificellement les NA par 0
    shortdata[is.na(shortdata)] <- 0
    if (redquali==TRUE){
      shortdata <- uniquecombs(shortdata)
      weightdata <- as.numeric(table(attr(shortdata,"index")))
    }
    colnames(shortdata) <- colnames(x)
    modalitynames <- list()
    for (j in 1:d){
      modalitynames[[j]] <- levels(x[,j])
      if (length(modalitynames[[j]]) != length(unique(x[which(is.na(x[,j])==FALSE),j])))
        stop(paste("The number of observed modalities is not equal to the number of levels for variable", colnames(x)[j]))
    }
    output <-  new("VSLCMdataCategorical", n=n, d=d, data=as.matrix(x), shortdata=shortdata, weightdata=weightdata, modalitynames=modalitynames)
  }else  if ( (length(idxcont) == d)){ 
    mat <- apply(x, 2, as.numeric)
    # construction des priors
    priors <- matrix(1, d, 4)
    priors[,4] <- 1/100
    priors[,3] <- colMeans(x, na.rm = T)
    colnames(priors) <- c("alpha", "beta", "lambda", "delta")
    ## Pour travailler avec Armadillo on rempli artificellement les NA par 0
    notNA <- (is.na(x)==FALSE)*1
    mat[is.na(mat)] <- 0
    colnames(mat) <-  colnames(x)
    colnames(notNA) <- colnames(x)
    output <-  new("VSLCMdataContinuous", n=n, d=d, data=mat, notNA=notNA, priors=priors)    
  }else  if ( (length(idxinte) == d)){ 
    mat <- apply(x, 2, as.numeric)
    # construction des priors
    priors <- matrix(1, d, 2)
    colnames(priors) <- c("alpha", "beta")
    ## Pour travailler avec Armadillo on rempli artificellement les NA par 0
    notNA <- (is.na(x)==FALSE)*1
    mat[is.na(mat)] <- 0
    colnames(mat) <-  colnames(x)
    colnames(notNA) <- colnames(x)
    output <-  new("VSLCMdataInteger", n=n, d=d, data=mat, notNA=notNA, priors=priors)    
  }else{
    output <- list(continuous=new("VSLCMdataContinuous"), integer=new("VSLCMdataInteger"), categorical=new("VSLCMdataCategorical"))
    if (length(idxcont) != 0){
      tmpdata <- data.frame(x[,idxcont])
      colnames(tmpdata) <- colnames(x)[idxcont]
      output$continuous <- VSLCMdata(tmpdata)
    }
    if (length(idxinte) != 0){
      tmpdata <- data.frame(x[,idxinte])
      colnames(tmpdata) <- colnames(x)[idxinte]
      output$integer <- VSLCMdata(tmpdata)
    }
    if (length(idxcat) != 0){
      tmpdata <- data.frame(x[,idxcat])
      colnames(tmpdata) <- colnames(x)[idxcat]      
      output$categorical <- VSLCMdata(tmpdata, redquali=FALSE)
    }
    
    output <- new("VSLCMdata", n=n, d=d, 
                  withContinuous=(length(idxcont) != 0),  withInteger=(length(idxinte) != 0), withCategorical=(length(idxcat) != 0),
                  dataContinuous=output$continuous, dataInteger=output$integer, dataCategorical=output$categorical,   var.names=colnames(x)
    )
  }
  return(output)
}
########################################################################################################################
## La fonction VSLCMdata permet de construire un objet de class S4 VSLCMdataContinuous ou VSLCMdataCategorical en fonction
## de la nature des variables
########################################################################################################################
VSLCMdataMixte <- function(x, redquali=TRUE){
  # Ajout d'un nom de variable si celui-ci est manquant
  if (is.null(colnames(x))) colnames(x) <- paste("X",1:ncol(x), sep="")
  if (max(table(colnames(x)))>1) stop("At least two variables have the same name!")
  n <- nrow(x)
  d <- ncol(x)
  # recherche des indices de variables numeric et factor
  type <- numeric()
  for (j in 1:d) type[j] <- class(x[,j])
  idxcont <- which(type=="numeric")
  idxinte <- which(type=="integer")
  idxcat <- which(type=="factor")
  if ((all(type %in% c("numeric", "integer", "factor"))==FALSE))
    stop("At least one variable is neither numeric, integer nor factor!")
  
  output <- list(continuous=new("VSLCMdataContinuous"), integer=new("VSLCMdataInteger"), categorical=new("VSLCMdataCategorical"))
  if (length(idxcont) != 0){
    tmpdata <- data.frame(x[,idxcont])
    colnames(tmpdata) <- colnames(x)[idxcont]
    output$continuous <- VSLCMdata(tmpdata)
  }
  if (length(idxinte) != 0){
    tmpdata <- data.frame(x[,idxinte])
    colnames(tmpdata) <- colnames(x)[idxinte]
    output$integer <- VSLCMdata(tmpdata)
  }
  if (length(idxcat) != 0){
    tmpdata <- data.frame(x[,idxcat])
    colnames(tmpdata) <- colnames(x)[idxcat]      
    output$categorical <- VSLCMdata(tmpdata, redquali=FALSE)
  }
  output <- new("VSLCMdata", n=n, d=d, 
                withContinuous=(length(idxcont) != 0),  withInteger=(length(idxinte) != 0), withCategorical=(length(idxcat) != 0),
                dataContinuous=output$continuous, dataInteger=output$integer, dataCategorical=output$categorical,   var.names=colnames(x)
  )
  return(output)
}

