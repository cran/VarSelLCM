
##'  Variable Selection for Model-Based Clustering of Mixed-Type Data Set with Missing Values
##'
##' Model-based clustering with variable selection and estimation of the number of clusters. Data to analyze can be continuous, categorical, integer or mixed. Moreover, missing values can occur and do not necessitate any pre-processing. Shiny application permits an easy interpretation of the results.
##'
##' \tabular{ll}{
##'   Package: \tab VarSelLCM\cr 
##'   Type: \tab Package\cr 
##'   Version: \tab 2.1.3\cr
##'   Date: \tab 2018-08-27\cr 
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
##' Standard tool methods (e.g., \link{summary}, \link{print}, \link{plot}, \link{coef}, \link{fitted}, \link{predict}...) are available for facilitating the interpretation.
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
##' Matthieu Marbac and Mohammed Sedki. Maintainer: Mohammed Sedki <mohammed.sedki@u-psud.fr>
##'
##' @references Marbac, M. and Sedki, M. (2017). Variable selection for model-based clustering using the integrated completed-data likelihood. Statistics and Computing, 27 (4), 1049-1063.
##' 
##' Marbac, M. and Patin, E. and Sedki, M. (2018). Variable selection for mixed data clustering: Application in human population genomics. Journal of classification, to appear.
##' 
##' @examples
##' \dontrun{
#' # Package loading
#' require(VarSelLCM)
#' 
#' # Data loading:
#' # x contains the observed variables
#' # z the known statu (i.e. 1: absence and 2: presence of heart disease)
#' data(heart)
#' ztrue <- heart[,"Class"]
#' x <- heart[,-13]
#' 
#' # Cluster analysis without variable selection
#' res_without <- VarSelCluster(x, 2, vbleSelec = FALSE, crit.varsel = "BIC")
#' 
#' # Cluster analysis with variable selection (with parallelisation)
#' res_with <- VarSelCluster(x, 2, nbcores = 2, initModel=40, crit.varsel = "BIC")
#' 
#' # Comparison of the BIC for both models:
#' # variable selection permits to improve the BIC
#' BIC(res_without)
#' BIC(res_with)
#' 
#' # Comparison of the partition accuracy. 
#' # ARI is computed between the true partition (ztrue) and its estimators
#' # ARI is an index between 0 (partitions are independent) and 1 (partitions are equals)
#' # variable selection permits to improve the ARI
#' # Note that ARI cannot be used for model selection in clustering, because there is no true partition
#' ARI(ztrue, fitted(res_without))
#' ARI(ztrue, fitted(res_with))
#' 
#' # Estimated partition
#' fitted(res_with)
#' 
#' # Estimated probabilities of classification
#' head(fitted(res_with, type="probability"))
#' 
#' # Summary of the probabilities of missclassification
#' plot(res_with, type="probs-class")
#' 
#' # Confusion matrices and ARI (only possible because the "true" partition is known).
#' # ARI is computed between the true partition (ztrue) and its estimators
#' # ARI is an index between 0 (partitions are independent) and 1 (partitions are equals)
#' # variable selection permits to improve the ARI
#' # Note that ARI cannot be used for model selection in clustering, because there is no true partition
#' # variable selection decreases the misclassification error rate
#' table(ztrue, fitted(res_without))
#' table(ztrue, fitted(res_with))
#' ARI(ztrue,  fitted(res_without))
#' ARI(ztrue, fitted(res_with))
#' 
#' # Summary of the best model
#' summary(res_with)
#' 
#' # Discriminative power of the variables (here, the most discriminative variable is MaxHeartRate)
#' plot(res_with)
#' 
#' # More detailed output
#' print(res_with)
#' 
#' # Print model parameter
#' coef(res_with)
#' 
#' # Boxplot for the continuous variable MaxHeartRate
#' plot(x=res_with, y="MaxHeartRate")
#' 
#' # Empirical and theoretical distributions of the most discriminative variable
#' # (to check that the distribution is well-fitted)
#' plot(res_with, y="MaxHeartRate", type="cdf")
#' 
#' # Summary of categorical variable
#' plot(res_with, y="Sex")
#' 
#' # Probabilities of classification for new observations 
#' predict(res_with, newdata = x[1:3,])
#' 
#' # Imputation by posterior mean for the first observation
#' not.imputed <- x[1,]
#' imputed <- VarSelImputation(res_with, x[1,], method = "sampling")
#' rbind(not.imputed, imputed)
#' 
#' # Opening Shiny application to easily see the results
#' VarSelShiny(res_with)
#' 
##' 
##' }
##' 
NULL

##' Statlog (Heart) Data Set
##' 
##'  This dataset is a heart disease database similar to a database already present in the repository (Heart Disease databases) but in a slightly different form.
##' 
##' 12 variables are used to cluster the observations 
##' \itemize{
#'  \item{age (integer)}
#'  \item{sex (binary)}
#'  \item{chest pain type (categorical with 4 levels)}
#'  \item{resting blood pressure (continuous)}
#'  \item{serum cholestoral in mg/dl (continuous)}
#'  \item{fasting blood sugar > 120 mg/dl (binary)}
#'  \item{resting electrocardiographic results (categorical with 3 levels)}
#'  \item{maximum heart rate achieved (continuous)}
#'  \item{exercise induced angina (binary)}
#'  \item{the slope of the peak exercise ST segment (categorical with 3 levels)}
#'  \item{number of major vessels  colored by flourosopy  (categorical with 4 levels)}
#'  \item{thal: 3 = normal; 6 = fixed defect; 7 = reversable defect (categorical with 3 levels)}
#'  }
#'  
#'  1 variable define a ''true'' partition: Absence (1) or presence (2) of heart disease 
##' 
##'
##' 
##' @references  UCI Machine Learning Repository [http://archive.ics.uci.edu/ml]. Irvine, CA: University of California, School of Information and Computer Science: http://archive.ics.uci.edu/ml/datasets/statlog+(heart)
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
##' @param x data.frame/matrix. Rows correspond to observations and columns correspond to variables. Continuous variables must be "numeric", count variables must be "integer" and categorical variables must be "factor"
##' @param gvals numeric. It defines number of components to consider.
##' @param vbleSelec logical. It indicates if a variable selection is done
##' @param crit.varsel character. It defines the information criterion used for model selection. Without variable selection, you can use one of the three criteria: "AIC", "BIC" and "ICL". With variable selection, you can use "AIC", BIC" and "MICL".
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
##' Marbac, M. and Patin, E. and Sedki, M. (2018). Variable selection for mixed data clustering: Application in human population genomics. Journal of Classification, to appear.
##' 
##' @examples
##' \dontrun{
#' # Package loading
#' require(VarSelLCM)
#' 
#' # Data loading:
#' # x contains the observed variables
#' # z the known statu (i.e. 1: absence and 2: presence of heart disease)
#' data(heart)
#' ztrue <- heart[,"Class"]
#' x <- heart[,-13]
#' 
#' # Cluster analysis without variable selection
#' res_without <- VarSelCluster(x, 2, vbleSelec = FALSE, crit.varsel = "BIC")
#' 
#' # Cluster analysis with variable selection (with parallelisation)
#' res_with <- VarSelCluster(x, 2, nbcores = 2, initModel=40, crit.varsel = "BIC")
#' 
#' # Comparison of the BIC for both models:
#' # variable selection permits to improve the BIC
#' BIC(res_without)
#' BIC(res_with)
#' 
#' # Confusion matrices and ARI (only possible because the "true" partition is known).
#' # ARI is computed between the true partition (ztrue) and its estimators
#' # ARI is an index between 0 (partitions are independent) and 1 (partitions are equals)
#' # variable selection permits to improve the ARI
#' # Note that ARI cannot be used for model selection in clustering, because there is no true partition
#' # variable selection decreases the misclassification error rate
#' table(ztrue, fitted(res_without))
#' table(ztrue, fitted(res_with))
#' ARI(ztrue,  fitted(res_without))
#' ARI(ztrue, fitted(res_with))
#'  
#' # Estimated partition
#' fitted(res_with)
#' 
#' # Estimated probabilities of classification
#' head(fitted(res_with, type="probability"))
#' 
#' # Summary of the probabilities of missclassification
#' plot(res_with, type="probs-class")
#' 
#' # Summary of the best model
#' summary(res_with)
#' 
#' # Discriminative power of the variables (here, the most discriminative variable is MaxHeartRate)
#' plot(res_with)
#' 
#' # More detailed output
#' print(res_with)
#' 
#' # Print model parameter
#' coef(res_with)
#' 
#' # Boxplot for the continuous variable MaxHeartRate
#' plot(x=res_with, y="MaxHeartRate")
#' 
#' # Empirical and theoretical distributions of the most discriminative variable 
#' # (to check that the distribution is well-fitted)
#' plot(res_with, y="MaxHeartRate", type="cdf")
#' 
#' # Summary of categorical variable
#' plot(res_with, y="Sex")
#' 
#' # Probabilities of classification for new observations 
#' predict(res_with, newdata = x[1:3,])
#' 
#' # Imputation by posterior mean for the first observation
#' not.imputed <- x[1,]
#' imputed <- VarSelImputation(res_with, x[1,], method = "sampling")
#' rbind(not.imputed, imputed)
#' 
#' # Opening Shiny application to easily see the results
#' VarSelShiny(res_with)
#' 
##' 
##' }
##' 
##' @export
##'
##'
VarSelCluster <- function(x, gvals, vbleSelec=TRUE, crit.varsel="BIC", initModel=50,  nbcores=1, discrim=rep(1,ncol(x)), nbSmall=250, iterSmall=20,  nbKeep=50, iterKeep=1000, tolKeep=10**(-6)){
  out <- list()
  for (g in 1:length(gvals))
    out[[g]] <- VarSelCluster.singleg(x, gvals[g], vbleSelec, crit.varsel, initModel,  nbcores, discrim, nbSmall, iterSmall,  nbKeep, iterKeep, tolKeep)
  if (crit.varsel=="BIC")  out <- out[[which.max(sapply(out, function(u) u@criteria@BIC))]]
  if (crit.varsel=="AIC") out <-   out[[which.max(sapply(out, function(u) u@criteria@AIC))]]
  if (crit.varsel=="ICL") out <-  out[[which.max(sapply(out, function(u) u@criteria@ICL))]]
  if (crit.varsel=="MICL") out <-  out[[which.max(sapply(out, function(u) u@criteria@MICL))]]
  out
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
setGeneric ( name= "MICLcomputation",  def = function(x, obj){ standardGeneric("MICLcomputation")})
## Pour les variables mixed
setMethod( f = "MICLcomputation", 
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

