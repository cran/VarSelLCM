

########################################################################################################################
## BIC extractor
########################################################################################################################
#'
#' BIC criterion.
#' 
#' This function gives the BIC criterion of an instance of \code{\linkS4class{VSLCMresults}}. 
#' BIC is computed according to the formula \deqn{BIC=log-likelihood - 0.5*\nu*log(n)} 
#' where  \eqn{\nu} denotes the number of parameters in the fitted model and \eqn{n} represents the sample size.  
#' 
#' @param object instance of \code{\linkS4class{VSLCMresults}}.
#' 
#' @name BIC
#' @rdname BIC-methods
#' @docType methods
#' @exportMethod BIC
#' @aliases BIC BIC,VSLCMresults-method
#' @references Schwarz, G. (1978). Estimating the dimension of a model. Annals of Statistics, 6(2), 461-464.
#' @examples
#' # Data loading:
#' data(heart)
#' 
#' # Cluster analysis without variable selection (number of clusters between 1 and 3)
#' res<- VarSelCluster(heart[,-13], 2, vbleSelec = FALSE)
#' 
#' # Get the BIC value
#' BIC(res)


setMethod(f="BIC",
          signature = c("VSLCMresults"),
          definition = function(object) object@criteria@BIC)

########################################################################################################################
## AIC extractor
########################################################################################################################
#'
#' AIC criterion.
#' 
#' This function gives the AIC criterion of an instance of \code{\linkS4class{VSLCMresults}}. 
#' AIC is computed according to the formula\deqn{AIC=log-likelihood - \nu} where  \eqn{\nu} denotes the number of parameters in the fitted model.
#' 
#' @param object instance of \code{\linkS4class{VSLCMresults}}.
#' 
#' @name AIC
#' @rdname AIC-methods
#' @docType methods
#' @exportMethod AIC
#' @aliases AIC AIC,VSLCMresults-method
#' @references Akaike, H. (1974), "A new look at the statistical model identification", IEEE Transactions on Automatic Control, 19 (6): 716-723.
#' @examples
#' # Data loading:
#' data(heart)
#' 
#' # Cluster analysis without variable selection
#' res <- VarSelCluster(heart[,-13], 2, vbleSelec = FALSE)
#' 
#' # Get the AIC value
#' AIC(res)

setMethod(f="AIC",
          signature = c("VSLCMresults"),
          definition = function(object) object@criteria@AIC)

########################################################################################################################
##  MICL extractor
########################################################################################################################
##' MICL criterion
##'
##' @description  
##' This function gives the MICL criterion for an instance of \code{\linkS4class{VSLCMresults}}.
##' 
##' @param object \code{\linkS4class{VSLCMresults}}
##' 
##' @references Marbac, M. and Sedki, M. (2017). Variable selection for model-based clustering using the integrated completed-data likelihood. Statistics and Computing, 27 (4), 1049-1063.
##' 
##' 
##' @examples
##' \dontrun{
##' # Data loading:
##' data("heart")
##' 
##' # Cluster analysis with variable selection
##' object <- VarSelCluster(heart[,-13], 2, vbleSelec = TRUE, crit.varsel = "MICL")
##' 
##' # Get the MICL value
##' MICL(object)
##' }
##' @export
##'
##'
MICL <- function(object){
  check.results(object)
  if (length(object@criteria@MICL)==0) stop("This criterion wasn't computed during model selection")
  object@criteria@MICL
}

########################################################################################################################
##  ICL extractor
########################################################################################################################
##' ICL criterion
##'
##' @description  
##' This function gives the ICL criterion for an instance of \code{\linkS4class{VSLCMresults}}.
##' 
##' @param object \code{\linkS4class{VSLCMresults}}
##' 
##' @references Biernacki, C., Celeux, G., and Govaert, G. (2000). Assessing a mixture model for clustering with the integrated completed likelihood. IEEE transactions on pattern analysis and machine intelligence, 22(7), 719-725.
##' 
##' 
##' @examples
##' # Data loading:
##' data(heart)
##' 
##' # Cluster analysis without variable selection
##' res <- VarSelCluster(heart[,-13], 2, vbleSelec = FALSE)
##' 
##' # Get the ICL value
##' ICL(res)
##' 
##' @export
##'
##'
ICL <- function(object){
  check.results(object)
  object@criteria@ICL
}

########################################################################################################################
## fitted 
########################################################################################################################
#'
#' Extract the partition or the probabilities of classification
#' 
#' @description  
#' This function returns the probabilities of classification or the partition among the observations of an instance of \code{\linkS4class{VSLCMresults}}.
#' 
#' @param object instance of \code{\linkS4class{VSLCMresults}}.
#' @param type the type of prediction: probability of classification (probability) or the partition (partition)
#'
#' 
#' @name fitted
#' @rdname fitted-methods
#' @docType methods
#' @exportMethod fitted
#' @aliases fitted fitted,VSLCMresults-method
#' @examples
#' # Data loading:
#' data(heart)
#' 
#' # Cluster analysis without variable selection (number of clusters between 1 and 3)
#' res <- VarSelCluster(heart[,-13], 2, vbleSelec = FALSE)
#' 
#' # Get the ICL value
#' fitted(res)

setMethod(f="fitted",
          signature = c("VSLCMresults"),
          definition = function(object, type="partition"){
            if (!(type %in% c("probability", "partition")))
              stop("type must be probability or partition")
            out <- object@partitions@zMAP
            if (type=="probability") out <- object@partitions@tik
            out
          }
)

########################################################################################################################
## fitted.values
########################################################################################################################
#'
#' Extract the partition or the probabilities of classification
#' 
#' @description  
#' This function returns the probabilities of classification or the partition among the observations of an instance of \code{\linkS4class{VSLCMresults}}.
#' 
#' @param object instance of \code{\linkS4class{VSLCMresults}}.
#' @param type the type of prediction: probability of classification (probability) or the partition (partition)
#'
#' 
#' @name fitted.values
#' @rdname fitted.values-methods
#' @docType methods
#' @exportMethod fitted.values
#' @aliases fitted.values fitted.values,VSLCMresults-method
#' @examples
#' # Data loading:
#' data(heart)
#' 
#' # Cluster analysis without variable selection (number of clusters between 1 and 3)
#' res <- VarSelCluster(heart[,-13], 2, vbleSelec = FALSE)
#' 
#' # Get the ICL value
#' fitted.values(res)

setMethod(f="fitted.values",
          signature = c("VSLCMresults"),
          definition = function(object, type="partition"){
            if (!(type %in% c("probability", "partition")))
              stop("type must be probability or partition")
            out <- object@partitions@zMAP
            if (type=="probability") out <- object@partitions@tik
            out
          }
          )


########################################################################################################################
## coef
########################################################################################################################
#'
#' Extract the parameters
#' 
#' @description  
#' This function returns an instance of class \code{\linkS4class{VSLCMparam}} which contains the model parameters.
#' 
#' @param object instance of \code{\linkS4class{VSLCMresults}}.
#'  
#' 
#' @name coef
#' @rdname coef-methods
#' @docType methods
#' @exportMethod coef
#' @aliases coef coef,VSLCMresults-method
#' 
#' @examples
#' # Data loading:
#' data(heart)
#' 
#' # Cluster analysis without variable selection (number of clusters between 1 and 3)
#' res  <- VarSelCluster(heart[,-13], 1:3, vbleSelec = FALSE)
#' 
#' # Get the ICL value
#' coef(res)

setMethod(f="coef",
          signature = c("VSLCMresults"),
          definition = function(object) object@param)

########################################################################################################################
## coefficients
########################################################################################################################
#'
#' Extract the parameters
#' 
#' @description  
#' This function returns an instance of class \code{\linkS4class{VSLCMparam}} which contains the model parameters.
#' 
#' @param object instance of \code{\linkS4class{VSLCMresults}}.
#' 
#' @name coefficients
#' @rdname coefficients-methods
#' @docType methods
#' @exportMethod coefficients
#' @aliases coefficients coefficients,VSLCMresults-method
#' @examples
#' # Data loading:
#' data(heart)
#' 
#' # Cluster analysis without variable selection (number of clusters between 1 and 3)
#' res  <- VarSelCluster(heart[,-13], 1:3, vbleSelec = FALSE)
#' 
#' # Get the ICL value
#' coefficients(res)

setMethod(f="coefficients",
          signature = c("VSLCMresults"),
          definition = function(object) object@param)
