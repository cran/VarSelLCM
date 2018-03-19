
##'  Variable Selection for Model-Based Clustering of Mixed-Type Data Set with Missing Values
##'
##' Model-based clustering with variable selection and estimation of the number of clusters. Data to analyze can be continuous, categorical, integer or mixed. Moreover, missing values can occur and do not necessitate any pre-processing. Shiny application permits an easy interpretation of the results.
##'
##' \tabular{ll}{
##'   Package: \tab VarSelLCM\cr 
##'   Type: \tab Package\cr 
##'   Version: \tab 2.1.0\cr
##'   Date: \tab 2018-03-14\cr 
##'   License: \tab GPL-3\cr  
##'   LazyLoad: \tab yes\cr
##'   URL:  \tab http://varsellcm.r-forge.r-project.org/\cr
##' }
##'
##' The main function to use is \link{VarSelCluster}. Function \link{VarSelCluster} carries out the model selection (according to AIC, BIC or MICL) and maximum likelihood estimation.
##' 
##' Function \link{VarSelShiny} runs a shiny application which permits an easy interpretation of the clustering results.
##' 
##' Function \link{VarSelImputation} permits the imputation of missing values by using the model parameters. 
##' 
##' Tool methods \link{summary}, \link{print} and \link{plot} are also available for facilitating the interpretation.
##' 
##' @name VarSelLCM-package
##' @aliases VarSelLCM
##' @rdname VarSelLCM-package
##' @docType package
##' @keywords package
##' @import parallel
##' @import Rcpp
##' @import methods
##' @import ggplot2
##' @import shiny
##' @importFrom mgcv uniquecombs
##' @importFrom graphics barplot mtext par
##' @importFrom stats dnorm dpois integrate sd runif pnorm ppois rnorm rpois
##' @useDynLib VarSelLCM
##'
##' @author
##' Matthieu Marbac and Mohammed Sedki Maintainer: Mohammed Sedki <mohammed.sedki@u-psud.fr>
##'
##' @references Marbac, M. and Sedki, M. (2017). Variable selection for model-based clustering using the integrated completed-data likelihood. Statistics and Computing, 27 (4), 1049-1063.
##' 
##' Marbac, M. and Patin, E. and Sedki, M. (2018). Variable selection for mixed data clustering: Application in human population genomics. Arxiv 1703.02293.
##' 
##' @examples
##' \dontrun{
##' # Package loading
##' require(VarSelLCM)
##' 
##' # Data loading:
##' # x contains the observed variables
##' # z the known statu (i.e. 1: absence and 2: presence of heart disease)
##' data(heart)
##' z <- heart[,"Class"]
##' x <- heart[,-13]
##' 
##' # Cluster analysis without variable selection
##' res_without <- VarSelCluster(x, 2, vbleSelec = FALSE)
##' 
##' # Cluster analysis with variable selection (with parallelisation)
##' res_with <- VarSelCluster(x, 2, nbcores = 2, initModel=40)
##' 
##' # Confusion matrices and ARI: variable selection decreases the misclassification error rate
##' print(table(z, res_without@partitions@zMAP))
##' print(table(z, res_with@partitions@zMAP))
##' ARI(z, res_without@partitions@zMAP)
##' ARI(z, res_with@partitions@zMAP)
##' 
##' # Summary of the best model
##' summary(res_with)
##' 
##' # Parameters of the best model
##' print(res_with)
##' 
##' # Opening Shiny application to easily see the results
##' VarSelShiny(res_with)
##' 
##' # Discriminative power of the variables (here, the most discriminative variable is MaxHeartRate)
##' plot(out, type="bar")
##' # Boxplot for continuous (or interger) variable
##' plot(out, y="MaxHeartRate", type="boxplot")
##'
##' # Empirical and theoretical distributions (to check that clustering is pertinent)
##' plot(out, y="MaxHeartRate", type="cdf")
##'
##'# Summary of categorical variable
##' plot(out, y="Sex")
##' 
##' # Summary of the probabilities of missclassification
##' plot(out, type="probs-class")
##' 
##' # Imputation by posterior mean for the first observation
##' not.imputed <- heart[1,-13]
##' imputed <- VarSelImputation(out)[1,]
##' rbind(not.imputed, imputed)
##' 
##' }
##' 
NULL

##' Statlog (Heart) Data Set
##' 
##'  This dataset is a heart disease database similar to a database already present in the repository (Heart Disease databases) but in a slightly different form.
##'  
##'  
##' 
##'
##' 
##' @references Website:http://archive.ics.uci.edu/ml/datasets/statlog+(heart)
##' @name heart
##' @docType data
##' @keywords datasets
##' 
##' @examples
##'   data(heart)
NULL

VarSelCluster.singleg <- function(x, g, vbleSelec, crit.varsel, initModel,  nbcores, discrim, nbSmall, iterSmall,  nbKeep, iterKeep, tolKeep){
  reference <- BuildS4Reference(x, g, initModel, vbleSelec, crit.varsel, TRUE, nbcores, discrim, nbSmall, iterSmall, nbKeep, iterKeep, tolKeep)
  if (g==1){
    reference <- withoutmixture(reference)
  }else{
    # Estimation du modele et/ou des parametres
    if (reference@strategy@parallel == FALSE)
      reference <- VarSelModelMLE(reference, 0)
    else{
      reference <- ParallelCriterion(reference, min(reference@strategy@initModel, nbcores))
      reference@strategy@vbleSelec <- vbleSelec
    }    
  }
  out <- DesignOutput(reference)
  out <- new("VSLCMresults", data=convertdata(out@data), criteria=out@criteria, partitions=out@partitions,
             model=out@model, strategy=out@strategy, param=convertparam(out@param))
  out@criteria@discrim <- sort(pvdiscrim(out), decreasing = TRUE)
  return(out)
}

###################################################################################
##' Variable selection and clustering.
##'
##' @description  
##' This function performs the model selection and the maximum likelihood estimation.
##' It can be used for clustering only (i.e., all the variables are assumed to be discriminative). In this case, you must specify the data to cluster (arg. x), the number of clusters (arg. g) and the option vbleSelec must be FALSE.
##' This function can also be used for variable selection in clustering. In this case, you must specify the data to analyse (arg. x), the number of clusters (arg. g) and the option vbleSelec must be TRUE. Variable selection can be done with BIC, MICL or AIC.
##'
##' @param x data.frame. Rows correspond to observations and columns correspond to variables. Continuous variables must be "numeric", count variables must be "integer" and categorical variables must be "factor"
##' @param gvals numeric. It defines number of components to consider.
##' @param vbleSelec logical. It indicates if a variable selection is done
##' @param crit.varsel character. It defines the information criterion used for the variable selection ("AIC", "BIC" or "MICL"; only used if vbleSelec=1)
##' @param initModel numeric. It gives the number of initializations of the alternated algorithm maximizing the MICL criterion (only used if crit.varsel="MICL")
##' @param nbcores numeric.  It defines the numerber of cores used by the alogrithm
##' @param discrim numeric. It indicates if each variable is discrimiative (1) or irrelevant (0) (only used if vbleSelec=0)
##' @param nbSmall numeric. It indicates  the number of SmallEM algorithms performed for the ML inference
##' @param iterSmall numeric. It indicates  the number of iterations for each SmallEM algorithm
##' @param nbKeep numeric. It indicates the number of chains used for the final EM algorithm
##' @param iterKeep numeric. It indicates the maximal number of iterations for each EM algorithm
##' @param tolKeep numeric. It indicates the maximal gap between two successive iterations of EM algorithm which stops the algorithm
##' 
##' @return Returns an instance of \linkS4class{VSLCMresults}.
##' 
##' @references Marbac, M. and Sedki, M. (2017). Variable selection for model-based clustering using the integrated completed-data likelihood. Statistics and Computing, 27 (4), 1049-1063.
##' 
##' Marbac, M. and Patin, E. and Sedki, M. (2018). Variable selection for mixed data clustering: Application in human population genomics. Arxiv 1703.02293.
##' 
##' @examples
##' \dontrun{
##' # Package loading
##' require(VarSelLCM)
##' 
##' # Data loading:
##' # x contains the observed variables
##' # z the known statu (i.e. 1: absence and 2: presence of heart disease)
##' data(heart)
##' z <- heart[,"Class"]
##' x <- heart[,-13]
##' 
##' # Cluster analysis without variable selection
##' res_without <- VarSelCluster(x, 2, vbleSelec = FALSE)
##' 
##' # Cluster analysis with variable selection (with parallelisation)
##' res_with <- VarSelCluster(x, 2, nbcores = 2, initModel=40)
##' 
##' # Confusion matrices and ARI: variable selection decreases the misclassification error rate
##' print(table(z, res_without@partitions@zMAP))
##' print(table(z, res_with@partitions@zMAP))
##' ARI(z, res_without@partitions@zMAP)
##' ARI(z, res_with@partitions@zMAP)
##' 
##' # Summary of the best model
##' summary(res_with)
##' 
##' # Opening Shiny application to easily see the results
##' VarSelShiny(res_with)
##' 
##' # Parameters of the best model
##' print(res_with)
##' 
##' # Discriminative power of the variables (here, the most discriminative variable is MaxHeartRate)
##' plot(out, type="bar")
##' # Boxplot for continuous (or interger) variable
##' plot(out, y="MaxHeartRate", type="boxplot")
##'
##' # Empirical and theoretical distributions (to check that clustering is pertinent)
##' plot(out, y="MaxHeartRate", type="cdf")
##'
##'# Summary of categorical variable
##' plot(out, y="Sex")
##' 
##' # Summary of the probabilities of missclassification
##' plot(out, type="probs-class")
##' 
##' # Imputation by posterior mean for the first observation
##' not.imputed <- heart[1,-13]
##' imputed <- VarSelImputation(out)[1,]
##' rbind(not.imputed, imputed)
##' 
##' }
##' @export
##'
##'
VarSelCluster <- function(x, gvals, vbleSelec=TRUE, crit.varsel="BIC", initModel=50,  nbcores=1, discrim=rep(1,ncol(x)), nbSmall=250, iterSmall=20,  nbKeep=50, iterKeep=1000, tolKeep=10**(-6)){
  out <- list()
  for (g in 1:length(gvals))
    out[[g]] <- VarSelCluster.singleg(x, gvals[g], vbleSelec, crit.varsel, initModel,  nbcores, discrim, nbSmall, iterSmall,  nbKeep, iterKeep, tolKeep)
  out[[which.max(sapply(out, function(u) u@criteria@BIC))]]
}

ParallelCriterion <- function(reference, nb.cpus){
  if (reference@strategy@vbleSelec==FALSE){
    reference@strategy@paramEstim <- TRUE
    reference@strategy@nbSmall <- ceiling(reference@strategy@nbSmall / nb.cpus)
    reference@strategy@nbKeep <- ceiling(reference@strategy@nbKeep / nb.cpus)
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(nb.cpus)
      common.objects <- c("reference","VarSelModelMLE", "OptimizeMICL")
      clusterEvalQ(cl, {require(VarSelLCM)})
      clusterExport(cl=cl, varlist = common.objects, envir = environment())
      reference <- parLapply(cl = cl, X  = as.list(rep(0, nb.cpus)), fun = function(g){VarSelModelMLE(reference,g)})
      stopCluster(cl)
    }else{
      reference <- mclapply(X = as.list(rep(0, nb.cpus)), FUN = VarSelModelMLE, obj=reference, mc.cores = nb.cpus)
    }
    # On conserve les parametres maximisant la vraisemblance
    tmploglike <- rep(NA, length(reference))
    for (it in 1:length(tmploglike)) {if (reference[[it]]@criteria@degeneracyrate!=1) tmploglike[it] <- reference[[it]]@criteria@loglikelihood}
    if (all(is.na(tmploglike))) tmploglike[1]=1
    reference <- reference[[which.max(tmploglike)]]
  }else{
    if (reference@strategy@crit.varsel=="MICL"){
      reference@strategy@paramEstim <- FALSE
      reference@strategy@initModel <- ceiling(reference@strategy@initModel / nb.cpus)       
      if(Sys.info()["sysname"] == "Windows"){
        cl <- makeCluster(nb.cpus)
        common.objects <- c("reference","VarSelModelMLE", "OptimizeMICL")
        clusterEvalQ(cl, {require(VarSelLCM)})
        clusterExport(cl=cl, varlist = common.objects, envir = environment())
        reference <- parLapply(cl = cl, X  = as.list(rep(0, nb.cpus)), fun = function(g){VarSelModelMLE(reference,g)})
        stopCluster(cl)          
      }else{
        reference <- mclapply(X = as.list(rep(0, nb.cpus)), FUN = VarSelModelMLE, obj=reference, mc.cores = nb.cpus)
      }
      # On conserve le meilleur modele au sens de MICL
      tmpMICL <- rep(NA, length(reference))
      for (it in 1:length(reference)) tmpMICL[it] <- reference[[it]]@criteria@MICL
      cvrate <- 0
      for (it in which(tmpMICL==max(tmpMICL))) cvrate <- cvrate + reference[[it]]@criteria@cvrate
      reference <- reference[[which.max(tmpMICL)]]
      reference@criteria@cvrate <- cvrate
      reference@strategy@paramEstim <- TRUE
      reference@strategy@vbleSelec <- FALSE
      reference <- ParallelCriterion(reference, nb.cpus)
    }else{
      reference@strategy@initModel <- ceiling(reference@strategy@initModel / nb.cpus)       
      if(Sys.info()["sysname"] == "Windows"){
        cl <- makeCluster(nb.cpus)
        common.objects <- c("reference","VarSelModelMLE", "OptimizeMICL")
        clusterEvalQ(cl, {require(VarSelLCM)})
        clusterExport(cl=cl, varlist = common.objects, envir = environment())
        reference <- parLapply(cl = cl, X  = as.list(rep(0, nb.cpus)), fun = function(g){VarSelModelMLE(reference,g)})
        stopCluster(cl)          
      }else{
        reference <- mclapply(X = as.list(rep(0, nb.cpus)), FUN = VarSelModelMLE, obj=reference, mc.cores = nb.cpus)
      }
      tmpIC <- rep(NA, length(reference))
      if (reference[[1]]@strategy@crit.varsel=="BIC") for (it in 1:length(reference)) tmpIC[it] <- reference[[it]]@criteria@BIC
      if (reference[[1]]@strategy@crit.varsel=="AIC") for (it in 1:length(reference)) tmpIC[it] <- reference[[it]]@criteria@AIC
      reference <- reference[[which.max(tmpIC)]]
    }
  }
  reference
}


BuildS4Reference <- function(x, g, initModel, vbleSelec, crit.varsel, paramEstim, nbcores, discrim, nbSmall, iterSmall, nbKeep, iterKeep, tolKeep){
  CheckInputs(x, g, initModel, vbleSelec, crit.varsel, discrim, paramEstim, nbcores, nbSmall, iterSmall, nbKeep, iterKeep, tolKeep)
  # Creation de l'objet S4 VSLCMstrategy contenant les parametres de reglage
  strategy <- VSLCMstrategy(initModel, nbcores, vbleSelec, crit.varsel, paramEstim, nbSmall, iterSmall, nbKeep, iterKeep, tolKeep)    
  # Creation de l'objet S4 VSLCMdataContinuous ou VSLCMdataCategorical
  if ((vbleSelec==FALSE) || (crit.varsel=="MICL"))
    data <- VSLCMdata(x)
  else
    data <- VSLCMdataMixte(x)  
  if (class(data) == "VSLCMdataContinuous")
    reference <- new("VSLCMresultsContinuous", data=data, criteria=InitCriteria(), model=new("VSLCMmodel",g=g, omega=discrim), strategy=strategy)
  else if (class(data) == "VSLCMdataInteger")
    reference <- new("VSLCMresultsInteger", data=data, criteria=InitCriteria(), model=new("VSLCMmodel",g=g, omega=discrim), strategy=strategy)
  else if (class(data) == "VSLCMdataCategorical")
    reference <- new("VSLCMresultsCategorical", data=data, criteria=InitCriteria(), model=new("VSLCMmodel",g=g, omega=discrim), strategy=strategy)
  else if (class(data) == "VSLCMdata")
    reference <- new("VSLCMresults", data=data, criteria=InitCriteria(), model=new("VSLCMmodel",g=g, omega=discrim), strategy=strategy)
  else
    stop("Problem in the data!")      
  return(reference)
}

# ParallelMICL <- function(reference, nb.cpus){
#   if (reference@strategy@crit.varsel == TRUE){
#       }
#   
#   return(reference)
# }
########################################################################################################################
## Fonctions principales du package, les seules accessibles par l'utilisateur sont VarSelCluster,
## Imputation (voir Imputation.R) et MICL
########################################################################################################################
setGeneric ( name= "MICL",  def = function(x, obj){ standardGeneric("MICL")})
## Pour les variables continues
setMethod( f = "MICL", 
           signature(x="data.frame", obj="VSLCMresultsContinuous"), 
           definition = function(x, obj){
             obj@strategy@crit.varsel <- TRUE
             obj@data  <- VSLCMdata(x)
             tmp <- ComputeMICL(obj, "Continuous")
             return(list(MICL=tmp@criteria@MICL, zOPT=tmp@partitions@zOPT+1))         
           }
)
## Pour les variables entieres
setMethod( f = "MICL", 
           signature(x="data.frame", obj="VSLCMresultsInteger"), 
           definition = function(x, obj){
             obj@strategy@crit.varsel <- TRUE
             # travail sur les donnees manquantes
             obj@data  <- VSLCMdata(x)
             tmp <- ComputeMICL(obj, "Integer")
             return(list(MICL=tmp@criteria@MICL, zOPT=tmp@partitions@zOPT+1))          
           }
)
## Pour les variables categorielles
setMethod( f = "MICL", 
           signature(x="data.frame", obj="VSLCMresultsCategorical"), 
           definition = function(x, obj){
             obj@strategy@crit.varsel <- TRUE
             obj@data  <- VSLCMdata(x)
             tmp <- ComputeMICL(obj, "Categorical")
             tmp@partitions@zOPT <-  1 + as.numeric(obj@partitions@zOPT[attr(obj@data@shortdata,"index")])
             return(list(MICL=tmp@criteria@MICL, zOPT=tmp@partitions@zOPT))       
           }
)
## Pour les variables mixed
setMethod( f = "MICL", 
           signature(x="data.frame", obj="VSLCMresults"), 
           definition = function(x, obj){
             obj@strategy@crit.varsel <- TRUE
             obj@data  <- VSLCMdata(x)
             tmp <- ComputeMICL(obj, "Mixed")
             return(list(MICL=tmp@criteria@MICL, zOPT=tmp@partitions@zOPT+1))     
           }
)


########################################################################################################################
## La fonction VarSelModelMLE permet d'effectuer l'estimation des parametres en considerant que les variables donnees
## dans le slot model de l'objet VSLCMresultsContinuous ou VSLCMresultsCategorical.
## Il appelle le code c++ et retourne un objet VSLCMresultsContinuous ou VSLCMresultsCategorical en fonction de la
## nature des donnees.
########################################################################################################################
setGeneric ( name= "VarSelModelMLE",  def = function(obj,it){ standardGeneric("VarSelModelMLE")})
## Pour les variables continues
setMethod( f = "VarSelModelMLE", 
           signature(obj="VSLCMresultsContinuous",it="numeric"), 
           definition = function(obj,it){
             if ((obj@strategy@vbleSelec==FALSE)||(obj@strategy@crit.varsel=="MICL")){
               reference <- OptimizeMICL(obj, "Continuous")               
             }else{
               stop("error")
             }
             return(reference)         
           }
)
## Pour les variables entiers
setMethod( f = "VarSelModelMLE", 
           signature(obj="VSLCMresultsInteger",it="numeric"), 
           definition = function(obj,it){
             if ((obj@strategy@vbleSelec==FALSE)||(obj@strategy@crit.varsel=="MICL")){
               reference <- OptimizeMICL(obj, "Integer")
             }else{
               stop("error")
             }
             return(reference)         
           }
)
## Pour les variables categorielles
setMethod( f = "VarSelModelMLE", 
           signature(obj="VSLCMresultsCategorical",it="numeric"), 
           definition = function(obj, it){
             if ((obj@strategy@vbleSelec==FALSE)||(obj@strategy@crit.varsel=="MICL")){
               reference <- OptimizeMICL(obj, "Categorical")
             }else{
               stop("error")
             }
             return(reference)           
           }
)
## Pour les variables mixed
setMethod( f = "VarSelModelMLE", 
           signature(obj="VSLCMresults",it="numeric"), 
           definition = function(obj, it){
             if ((obj@strategy@vbleSelec==FALSE)||(obj@strategy@crit.varsel=="MICL")){
               reference <- OptimizeMICL(obj, "Mixed")
             }else{
               pen <- 0
               if (obj@strategy@crit.varsel=="AIC") pen <- 1
               if (obj@strategy@crit.varsel=="BIC") pen <- 0.5*log(obj@data@n)
               reference <- OptimizePenLike(obj, pen)
             }
             return(reference)           
           }
)

