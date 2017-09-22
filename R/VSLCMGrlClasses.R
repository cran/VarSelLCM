########################################################################################################################
## Differentes classes S4 accessibles a l'utilisateur et etant des slots des classes S4
## VSLCMresultsContinuous et/ou VSLCMresultsCategorical
########################################################################################################################

########################################################################################################################
## Classe S4 VSLCMcriteria contenant la logvraisemblance (loglikelihood), la valeur des criteres BIC, ICL et MICL
########################################################################################################################
###################################################################################
##' Constructor of [\code{\linkS4class{VSLCMcriteria}}] class
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
##' }
##'
##' @examples
##'   getSlots("VSLCMcriteria")
##'
##' @name VSLCMcriteria-class
##' @rdname VSLCMcriteria-class
##' @exportClass VSLCMcriteria
setClass(Class = "VSLCMcriteria", 
         representation = representation(loglikelihood="numeric", AIC="numeric", BIC="numeric", ICL="numeric", MICL="numeric", nbparam="numeric", cvrate="numeric", degeneracyrate="numeric"), 
         prototype = prototype(loglikelihood=numeric(), AIC=numeric(), BIC=numeric(), ICL=numeric(), MICL=numeric(), nbparam=numeric(), cvrate=numeric(), degeneracyrate=numeric() )
)

InitCriteria <- function()
  new("VSLCMcriteria", loglikelihood=-Inf, AIC=-Inf, BIC=-Inf, ICL=-Inf, MICL=-Inf)

########################################################################################################################
## Classe S4 VSLCMpartitions contenant la partition MAP (zMAP), la partition zstar (zOPT) et la partition floue (tik)
########################################################################################################################
###################################################################################
##' Constructor of [\code{\linkS4class{VSLCMpartitions}}] class
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
##' Constructor of [\code{\linkS4class{VSLCMstrategy}}] class
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
##' Constructor of [\code{\linkS4class{VSLCMmodel}}] class
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

########################################################################################################################
## Classe S4 VSLCMparamContinuous contenant les proportions (pi), les moyennes (mu) et les ecrat-types (sd)
########################################################################################################################
###################################################################################
##' Constructor of [\code{\linkS4class{VSLCMparamContinuous}}] class
##'
##'  
##' \describe{
##'   \item{pi}{numeric. Proportions of the mixture components.}
##'   \item{mu}{matrix. Mean for each component (column) and each variable (row).}
##'   \item{sd}{matrix. Standard deviation for each component (column) and each variable (row).}
##' }
##'
##' @examples
##'   getSlots("VSLCMparamContinuous")
##'
##' @name VSLCMparamContinuous-class
##' @rdname VSLCMparamContinuous-class
##' @exportClass VSLCMparamContinuous
setClass(
  Class = "VSLCMparamContinuous", 
  representation = representation(pi="numeric", mu="matrix", sd="matrix"), 
  prototype = prototype(pi=numeric(), mu=matrix(), sd=matrix())
)

########################################################################################################################
## Classe S4 VSLCMparamInteger contenant les proportions (pi), et les parameters (lambda)
########################################################################################################################
###################################################################################
##' Constructor of [\code{\linkS4class{VSLCMparamInteger}}] class
##'
##'  
##' \describe{
##'   \item{pi}{numeric. Proportions of the mixture components.}
##'   \item{lambda}{matrix. Mean for each component (column) and each variable (row).}
##' }
##'
##' @examples
##'   getSlots("VSLCMparamInteger")
##'
##' @name VSLCMparamInteger-class
##' @rdname VSLCMparamInteger-class
##' @exportClass VSLCMparamInteger
setClass(
  Class = "VSLCMparamInteger", 
  representation = representation(pi="numeric", lambda="matrix"), 
  prototype = prototype(pi=numeric(), lambda=matrix())
)

########################################################################################################################
## Classe S4 VSLCMparamCategorical contenant les proportions (pi), les probas (alpha)
########################################################################################################################
###################################################################################
##' Constructor of [\code{\linkS4class{VSLCMparamCategorical}}] class
##'
##'  
##' \describe{
##'   \item{pi}{numeric. Proportions of the mixture components.}
##'   \item{alpha}{list. Parameters of the multinomial distributions.}
##' }
##'
##' @examples
##'   getSlots("VSLCMparamCategorical")
##'
##' @name VSLCMparamCategorical-class
##' @rdname VSLCMparamCategorical-class
##' @exportClass VSLCMparamCategorical
setClass(
  Class = "VSLCMparamCategorical", 
  representation = representation(pi="numeric", alpha="list"), 
  prototype = prototype(pi=numeric(), alpha=list())
)


########################################################################################################################
## Classe S4 VSLCMparamMixed contenant les parametres continus et categoriels
########################################################################################################################
###################################################################################
##' Constructor of [\code{\linkS4class{VSLCMparamMixed}}] class
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
##'   getSlots("VSLCMparamMixed")
##'
##' @name VSLCMparamMixed-class
##' @rdname VSLCMparamMixed-class
##' @exportClass VSLCMparamMixed
setClass(
  Class = "VSLCMparamMixed", 
  representation = representation(pi="numeric", paramContinuous="VSLCMparamContinuous", paramInteger="VSLCMparamInteger", paramCategorical="VSLCMparamCategorical"), 
  prototype = prototype(pi=numeric(), paramContinuous=new("VSLCMparamContinuous"), paramInteger=new("VSLCMparamInteger"), paramCategorical=new("VSLCMparamCategorical"))
)

########################################################################################################################
## Classe S4 VSLCMresultsContinuous
########################################################################################################################.
###################################################################################
##' Constructor of [\code{\linkS4class{VSLCMresultsContinuous}}] class
##'
##'  
##' \describe{
##'   \item{data}{\linkS4class{VSLCMdataContinuous}. Results relied to the data.}
##'   \item{criteria}{\linkS4class{VSLCMcriteria}. Results relied to the information criteria.}
##'   \item{partitions}{\linkS4class{VSLCMpartitions}. Results relied to the partitions.}
##'   \item{model}{\linkS4class{VSLCMmodel}. Results relied to the selected model.}
##'   \item{strategy}{\linkS4class{VSLCMstrategy}. Results relied to the tune parameters.}
##'   \item{param}{\linkS4class{VSLCMparamContinuous}. Results relied to the parameters.}
##' }
##'
##' @examples
##'   getSlots("VSLCMresultsContinuous")
##'
##' @name VSLCMresultsContinuous-class
##' @rdname VSLCMresultsContinuous-class
##' @exportClass VSLCMresultsContinuous
setClass(
  Class = "VSLCMresultsContinuous", 
  representation = representation(data="VSLCMdataContinuous", criteria="VSLCMcriteria", partitions="VSLCMpartitions",
                                  model="VSLCMmodel", strategy="VSLCMstrategy", param="VSLCMparamContinuous"), 
  prototype = prototype(data=new("VSLCMdataContinuous"), criteria=new("VSLCMcriteria"), partitions=new("VSLCMpartitions"),
                        model=new("VSLCMmodel"), strategy=new("VSLCMstrategy"), param=new("VSLCMparamContinuous"))
)


########################################################################################################################
## Classe S4 VSLCMresultsInteger
########################################################################################################################
###################################################################################
##' Constructor of [\code{\linkS4class{VSLCMresultsInteger}}] class
##'
##'  
##' \describe{
##'   \item{data}{\linkS4class{VSLCMdataInteger}. Results relied to the data.}
##'   \item{criteria}{\linkS4class{VSLCMcriteria}. Results relied to the information criteria.}
##'   \item{partitions}{\linkS4class{VSLCMpartitions}. Results relied to the partitions.}
##'   \item{model}{\linkS4class{VSLCMmodel}. Results relied to the selected model.}
##'   \item{strategy}{\linkS4class{VSLCMstrategy}. Results relied to the tune parameters.}
##'   \item{param}{\linkS4class{VSLCMparamInteger}. Results relied to the parameters.}
##' }
##'
##' @examples
##'   getSlots("VSLCMresultsInteger")
##'
##' @name VSLCMresultsInteger-class
##' @rdname VSLCMresultsInteger-class
##' @exportClass VSLCMresultsInteger
setClass(
  Class = "VSLCMresultsInteger", 
  representation = representation(data="VSLCMdataInteger", criteria="VSLCMcriteria", partitions="VSLCMpartitions",
                                  model="VSLCMmodel", strategy="VSLCMstrategy", param="VSLCMparamInteger"), 
  prototype = prototype(data=new("VSLCMdataInteger"), criteria=new("VSLCMcriteria"), partitions=new("VSLCMpartitions"),
                        model=new("VSLCMmodel"), strategy=new("VSLCMstrategy"), param=new("VSLCMparamInteger"))
)

########################################################################################################################
## Classe S4 VSLCMresultsCategorical
########################################################################################################################
###################################################################################
##' Constructor of [\code{\linkS4class{VSLCMresultsCategorical}}] class
##'
##'  
##' \describe{
##'   \item{data}{\linkS4class{VSLCMdataCategorical}. Results relied to the data.}
##'   \item{criteria}{\linkS4class{VSLCMcriteria}. Results relied to the information criteria.}
##'   \item{partitions}{\linkS4class{VSLCMpartitions}. Results relied to the partitions.}
##'   \item{model}{\linkS4class{VSLCMmodel}. Results relied to the selected model.}
##'   \item{strategy}{\linkS4class{VSLCMstrategy}. Results relied to the tune parameters.}
##'   \item{param}{\linkS4class{VSLCMparamCategorical}. Results relied to the parameters.}
##' }
##'
##' @examples
##'   getSlots("VSLCMresultsCategorical")
##'
##' @name VSLCMresultsCategorical-class
##' @rdname VSLCMresultsCategorical-class
##' @exportClass VSLCMresultsCategorical
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
##' Constructor of [\code{\linkS4class{VSLCMresultsMixed}}] class
##'
##'  
##' \describe{
##'   \item{data}{\linkS4class{VSLCMdataMixed}. Results relied to the data.}
##'   \item{criteria}{\linkS4class{VSLCMcriteria}. Results relied to the information criteria.}
##'   \item{partitions}{\linkS4class{VSLCMpartitions}. Results relied to the partitions.}
##'   \item{model}{\linkS4class{VSLCMmodel}. Results relied to the selected model.}
##'   \item{strategy}{\linkS4class{VSLCMstrategy}. Results relied to the tune parameters.}
##'   \item{param}{\linkS4class{VSLCMparamMixed}. Results relied to the parameters.}
##' }
##'
##' @examples
##'   getSlots("VSLCMresultsMixed")
##'
##' @name VSLCMresultsMixed-class
##' @rdname VSLCMresultsMixed-class
##' @exportClass VSLCMresultsMixed
setClass(
  Class = "VSLCMresultsMixed", 
  representation = representation(data="VSLCMdataMixed", criteria="VSLCMcriteria", partitions="VSLCMpartitions",
                                  model="VSLCMmodel", strategy="VSLCMstrategy", param="VSLCMparamMixed"), 
  prototype = prototype(data=new("VSLCMdataMixed"), criteria=new("VSLCMcriteria"), partitions=new("VSLCMpartitions"),
                        model=new("VSLCMmodel"), strategy=new("VSLCMstrategy"), param=new("VSLCMparamMixed"))
)