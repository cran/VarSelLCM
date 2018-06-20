
########################################################################################################################
## Differentes classes S4 accessibles a l'utilisateur et etant des slots des classes S4
## VSLCMresultsContinuous et/ou VSLCMresultsCategorical
########################################################################################################################

########################################################################################################################
## Classe S4 VSLCMcriteria contenant la logvraisemblance (loglikelihood), la valeur des criteres BIC, ICL et MICL
########################################################################################################################
###################################################################################
##' Constructor of \code{\linkS4class{VSLCMcriteria}} class
##'
##'  
##' \describe{
##'   \item{loglikelihood}{numeric. Log-likelihood}
##'   \item{AIC}{numeric. Value of the AIC criterion.}
##'   \item{BIC}{numeric. Value of the BIC criterion.}
##'   \item{ICL}{numeric. Value of the ICL criterion.}
##'   \item{MICL}{numeric. Value of the MICL criterion.}
##'   \item{nbparam}{integer. Number of parameters.}
##'   \item{cvrate}{numeric.  Rate of convergence of the alternated algorithm for optimizing the MICL criterion.}
##'   \item{degeneracyrate}{numeric.  Rate of degeneracy for the selected model.}
##'   \item{discrim}{numeric.  Discriminative power of each variable.}
##' }
##'
##' @examples
##'   getSlots("VSLCMcriteria")
##'
##' @name VSLCMcriteria-class
##' @rdname VSLCMcriteria-class
##' @exportClass VSLCMcriteria
setClass(Class = "VSLCMcriteria", 
         representation = representation(loglikelihood="numeric", AIC="numeric", BIC="numeric", ICL="numeric", MICL="numeric", nbparam="numeric", cvrate="numeric", degeneracyrate="numeric", discrim="numeric"), 
         prototype = prototype(loglikelihood=numeric(), AIC=numeric(), BIC=numeric(), ICL=numeric(), MICL=numeric(), nbparam=numeric(), cvrate=numeric(), degeneracyrate=numeric(), discrim=numeric())
)

InitCriteria <- function()
  new("VSLCMcriteria", loglikelihood=-Inf, AIC=-Inf, BIC=-Inf, ICL=-Inf, MICL=-Inf)

########################################################################################################################
## Classe S4 VSLCMpartitions contenant la partition MAP (zMAP), la partition zstar (zOPT) et la partition floue (tik)
########################################################################################################################
###################################################################################
##' Constructor of \code{\linkS4class{VSLCMpartitions}} class
##'
##'  
##' \describe{
##'   \item{zMAP}{numeric. A vector indicating the class membership of each individual by using the MAP rule computed for the best model with its maximum likelihood estimates.}
##'   \item{zOPT}{numeric. Partition maximizing the integrated complete-data likelihood of the selected model.}
##'   \item{tik}{numeric. Fuzzy partition computed for the best model with its maximum likelihood estimates.}
##' }
##'
##' @examples
##'   getSlots("VSLCMpartitions")
##'
##' @name VSLCMpartitions-class
##' @rdname VSLCMpartitions-class
##' @exportClass VSLCMpartitions
setClass(
  Class = "VSLCMpartitions", 
  representation = representation(zMAP="numeric", zOPT="numeric", tik="matrix"), 
  prototype = prototype(zMAP=numeric(), zOPT=numeric(), tik=matrix(0,0,0))
)

########################################################################################################################
## Classe S4 VSLCMstrategy contenant les parametres de reglages detailles dans VarSELLCMmixte.R
########################################################################################################################
###################################################################################
##' Constructor of \code{\linkS4class{VSLCMstrategy}} class
##'
##'  
##' \describe{
##'   \item{initModel}{numeric. Number of initialisations for the model selection algorithm.}
##'   \item{vbleSelec}{logical. It indicates if the selection of the variables is performed.}
##'   \item{paramEstim}{logical. It indicates if the parameter estimation is performed.}
##'   \item{parallel}{logical. It indicates  if a parallelisation is done.}
##'   \item{nbSmall}{numeric. It indicates the number of small EM.}
##'   \item{iterSmall}{numeric. It indicates the number of iteration for the small EM}
##'   \item{nbKeep}{numeric. It indicates the number of chains kept for the EM.}
##'   \item{iterKeep}{numeric. It indicates the maximum number of iteration for the EM.}
##'   \item{tolKeep}{numeric. It indicates the value of the difference between successive iterations of EM stopping the EM.}
##' }
##'
##' @examples
##'   getSlots("VSLCMstrategy")
##'
##' @name VSLCMstrategy-class
##' @rdname VSLCMstrategy-class
##' @exportClass VSLCMstrategy
##'
setClass(
  Class = "VSLCMstrategy", 
  representation = representation(initModel="numeric", vbleSelec="logical", crit.varsel="character", paramEstim="logical", parallel="logical",
    nbSmall="numeric", iterSmall="numeric", nbKeep="numeric", iterKeep="numeric", tolKeep="numeric"), 
  prototype = prototype(initModel=numeric(),  vbleSelec=logical(), crit.varsel=character(), paramEstim=logical(), parallel=logical(),
    nbSmall=numeric(), iterSmall=numeric(), nbKeep=numeric(), iterKeep=numeric(), tolKeep=numeric())
) 

## Constructeur de la classe S4 VSLCMstrategy
VSLCMstrategy <- function(initModel, nbcores, vbleSelec, crit.varsel, paramEstim, nbSmall, iterSmall, nbKeep, iterKeep, tolKeep){
  new("VSLCMstrategy",
      initModel=initModel,
      parallel=(nbcores>1),
      vbleSelec=vbleSelec,
      crit.varsel=crit.varsel,
      paramEstim=paramEstim, 
      nbSmall=nbSmall, 
      iterSmall=iterSmall, 
      nbKeep=min(nbKeep, nbSmall),
      iterKeep=iterKeep, 
      tolKeep=tolKeep)
}


########################################################################################################################
## Classe S4 VSLCMmodel contenant le nombre de classes (g) et le role des variables (omega)
########################################################################################################################
###################################################################################
##' Constructor of \code{\linkS4class{VSLCMmodel}} class
##'
##'  
##' \describe{
##'   \item{g}{numeric. Number of components.}
##'   \item{omega}{logical. Vector indicating if each variable is irrelevant (1) or not (0) to the clustering.}
##'   \item{names.relevant}{character. Names of the relevant variables.}
##'   \item{names.irrelevant}{character. Names of the irrelevant variables.}
##' }
##'
##' @examples
##'   getSlots("VSLCMmodel")
##'
##' @name VSLCMmodel-class
##' @rdname VSLCMmodel-class
##' @exportClass VSLCMmodel
##'
setClass(
  Class = "VSLCMmodel", 
  representation = representation(g="numeric", omega="numeric", names.relevant="character", names.irrelevant="character"), 
  prototype = prototype(g=numeric(), omega=numeric(), names.relevant=character(), names.irrelevant=character())
)

check.results <- function(obj){
  if (class(obj)!="VSLCMresults") stop("Results must be an instance of VSLCMresults returned by the function VarSelCluster of R package VarSelLCM")
}