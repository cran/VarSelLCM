########################################################################################################################
## La fonction DesignOutput permet de mettre en forme les parametres en fonction de leur nature VSLCMresultsContinuous
## ou VSLCMresultsCategorical. Elle est appelee a la fin de l'estimation des parametres
########################################################################################################################
setGeneric ( name= "DesignOutput",  def = function(reference){ standardGeneric("DesignOutput")})

## Cas de variables continues
setMethod( f = "DesignOutput", 
           signature(reference="VSLCMresultsContinuous"), 
           definition = function(reference){
             reference@model@omega <-  as.numeric(reference@model@omega)
             names(reference@model@omega) <- colnames(reference@data@data)
             if (reference@model@g>1){
               if (reference@strategy@vbleSelec==FALSE){
                 reference@partitions@zOPT <- numeric()
                 reference@criteria@MICL <- numeric()
               }
               # On remet les valeurs manquantes
               for (j in 1:reference@data@d){
                 if (any(reference@data@notNA[,j]==0))
                   reference@data@data[which(reference@data@notNA[,j]==0),j] <- NA
               }
               if (reference@strategy@paramEstim){
                 if (reference@criteria@degeneracyrate != 1){
                   rownames(reference@param@mu)  <-   colnames(reference@data@data)
                   rownames(reference@param@sd)  <-   colnames(reference@data@data)
                   if (reference@criteria@degeneracyrate>0.5)
                     warning(paste("The rate of degeneracy for the EM algorithm is", reference@criteria@degeneracyrate ), call. = FALSE)
                   reference@param@pi <- as.numeric(reference@param@pi)
                   names(reference@param@pi) <-  paste("class-",1:length(reference@param@pi),sep="")
                   colnames(reference@param@mu) <-  paste("class-",1:length(reference@param@pi),sep="")
                   colnames(reference@param@sd) <-  paste("class-",1:length(reference@param@pi),sep="")
                   
                   reference@partitions@zMAP <- as.numeric(reference@partitions@zMAP) + 1
                   reference@partitions@zOPT <- as.numeric(reference@partitions@zOPT) + 1
                   colnames(reference@partitions@tik) <-  paste("class-",1:reference@model@g,sep="")
                   nbparam <- (reference@model@g-1 + reference@model@g*2*sum(reference@model@omega) + 2*sum(1-reference@model@omega))
                   reference@criteria@nbparam <- nbparam
                   reference@criteria@BIC <- reference@criteria@loglikelihood  - 0.5*nbparam*log(reference@data@n)
                   reference@criteria@AIC <- reference@criteria@loglikelihood - nbparam
                   reference@criteria@ICL <- ICLcontinuous(reference) 
                 }else{
                   warning("All the models got error (degeneracy)", call. = FALSE)
                 }
               }
             }
             names(reference@criteria@nbparam) <- NULL  
             if (any(reference@model@omega==1))
               reference@model@names.relevant <- as.character(names(reference@model@omega)[which(reference@model@omega==1)])
             if (any(reference@model@omega==0))
               reference@model@names.irrelevant <- names(reference@model@omega)[which(reference@model@omega==0)]
             
             return(reference)
           }
)

## Cas de variables entieres
setMethod( f = "DesignOutput", 
           signature(reference="VSLCMresultsInteger"), 
           definition = function(reference){
             reference@model@omega <-  as.numeric(reference@model@omega)
             names(reference@model@omega) <- colnames(reference@data@data)
             if (reference@model@g>1){
               if (reference@strategy@vbleSelec==FALSE){
                 reference@partitions@zOPT <- numeric()
                 reference@criteria@MICL <- numeric()
               }
               # On remet les valeurs manquantes
               for (j in 1:reference@data@d){
                 if (any(reference@data@notNA[,j]==0))
                   reference@data@data[which(reference@data@notNA[,j]==0),j] <- NA
               }
               if (reference@strategy@paramEstim){
                 if (reference@criteria@degeneracyrate != 1){
                   rownames(reference@param@lambda)  <-   colnames(reference@data@data)
                   if (reference@criteria@degeneracyrate>0.5)
                     warning(paste("The rate of degeneracy for the EM algorithm is", reference@criteria@degeneracyrate ), call. = FALSE)
                   reference@param@pi <- as.numeric(reference@param@pi)
                   names(reference@param@pi) <-  paste("class-",1:length(reference@param@pi),sep="")
                   colnames(reference@param@lambda) <-  paste("class-",1:length(reference@param@pi),sep="")           
                   reference@partitions@zMAP <- as.numeric(reference@partitions@zMAP) + 1
                   reference@partitions@zOPT <- as.numeric(reference@partitions@zOPT) + 1
                   colnames(reference@partitions@tik) <-  paste("class-",1:reference@model@g,sep="")
                   nbparam <- (reference@model@g-1 + reference@model@g*sum(reference@model@omega) + sum(1-reference@model@omega))
                   reference@criteria@nbparam <- nbparam
                   reference@criteria@BIC <- reference@criteria@loglikelihood  - 0.5*nbparam*log(reference@data@n)
                   reference@criteria@ICL <- ICLinteger(reference) 
                   reference@criteria@AIC <- reference@criteria@loglikelihood - nbparam
                 }else{
                   warning("All the models got error (degeneracy)", call. = FALSE)
                 }
               }
             }
             names(reference@criteria@nbparam) <- NULL
             if (any(reference@model@omega==1))
               reference@model@names.relevant <- as.character(names(reference@model@omega)[which(reference@model@omega==1)])
             if (any(reference@model@omega==0))
               reference@model@names.irrelevant <- names(reference@model@omega)[which(reference@model@omega==0)]
             
             return(reference)
           }
)
## Cas des variables categorielles
setMethod( f = "DesignOutput", 
           signature(reference="VSLCMresultsCategorical"), 
           definition = function(reference){
             reference@model@omega <-  as.numeric(reference@model@omega)
             names(reference@model@omega) <- colnames(reference@data@data)
             
             if (reference@model@g>1){
               if (reference@strategy@paramEstim == TRUE){
                 if (reference@strategy@vbleSelec==FALSE){
                   reference@partitions@zOPT <- numeric()
                   reference@criteria@MICL <- numeric()
                 }else{
                   reference@partitions@zOPT <- 1 + as.numeric(reference@partitions@zOPT[attr(reference@data@shortdata,"index")])
                 }
                 reference@partitions@zMAP <- 1 + as.numeric(reference@partitions@zMAP[attr(reference@data@shortdata,"index")])
                 # Attention zOPT correspond aux profiles, on repasse donc au niveau des individus
                 reference@param@pi <- as.numeric(reference@param@pi)
                 names(reference@param@pi) <- paste("class-",1:length(reference@param@pi),sep="")
                 for (j in 1:reference@data@d){
                   reference@param@alpha[[j]] <- matrix(reference@param@alpha[[j]], nrow = reference@model@g)
                   rownames(reference@param@alpha[[j]]) <- paste("class-",1:length(reference@param@pi),sep="")
                   colnames(reference@param@alpha[[j]]) <- reference@data@modalitynames[[j]]
                 }
                 names(reference@param@alpha) <- colnames(reference@data@shortdata)
                 if (reference@model@g>1)
                   reference@partitions@tik <- reference@partitions@tik[attr(reference@data@shortdata,"index"),] 
                 colnames(reference@partitions@tik)=paste("class-",1:reference@model@g,sep="")
                 nbparam <- (reference@model@g-1)
                 for (j in 1:reference@data@d)
                   nbparam <- nbparam + (length(reference@data@modalitynames[[j]])-1)*(1 + (reference@model@g-1)*reference@model@omega[j])
                 reference@criteria@nbparam <- nbparam
                 reference@criteria@BIC <- reference@criteria@loglikelihood  - 0.5*nbparam*log(reference@data@n)
                 reference@criteria@ICL <- ICLcategorical(reference) 
                 reference@criteria@AIC <- reference@criteria@loglikelihood - nbparam
               }
               
               # On remet les valeurs manquantes
               if (any(reference@data@shortdata == 0)){
                 for (j in 1:ncol(reference@data@shortdata))
                   reference@data@shortdata[which(reference@data@shortdata[,j] == 0), ] <- NA
               }
             }
             names(reference@criteria@nbparam) <- NULL
             if (any(reference@model@omega==1))
               reference@model@names.relevant <- as.character(names(reference@model@omega)[which(reference@model@omega==1)])
             if (any(reference@model@omega==0))
               reference@model@names.irrelevant <- names(reference@model@omega)[which(reference@model@omega==0)]
             
             return(reference)
           }
)

## Cas des variables mixed
setMethod( f = "DesignOutput", 
           signature(reference="VSLCMresultsMixed"), 
           definition = function(reference){
             reference@model@omega <-  as.numeric(reference@model@omega)
             namestmp <- numeric()
             if (reference@data@withContinuous) namestmp <- c(namestmp, colnames(reference@data@dataContinuous@data))
             if (reference@data@withInteger) namestmp <- c(namestmp, colnames(reference@data@dataInteger@data))
             if (reference@data@withCategorical) namestmp <- c(namestmp, colnames(reference@data@dataCategorical@data))
             names(reference@model@omega) <- namestmp
             if (reference@model@g>1){
               # On remet les valeurs manquantes
                              
               if (reference@strategy@paramEstim == TRUE){
                 if (reference@strategy@vbleSelec==FALSE){
                   reference@partitions@zOPT <- numeric()
                   reference@criteria@MICL <- numeric()
                   reference@partitions@zMAP <- as.numeric(reference@partitions@zMAP) + 1
                 }else{
                   reference@partitions@zMAP <- as.numeric(reference@partitions@zMAP) + 1
                   reference@partitions@zOPT <- as.numeric(reference@partitions@zOPT) + 1
                 }
                 
                 if (reference@criteria@degeneracyrate != 1){
                   if (reference@criteria@degeneracyrate>0.5)
                     warning(paste("The rate of degeneracy for the EM algorithm is", reference@criteria@degeneracyrate ), call. = FALSE)
                   
                   reference@param@pi <- as.numeric(reference@param@pi)
                   names(reference@param@pi) <- paste("class-",1:length(reference@param@pi),sep="")
                   reference@partitions@tik <- reference@partitions@tik 
                   colnames(reference@partitions@tik ) <- paste("class-",1:reference@model@g,sep="")
                   nbparam <- reference@model@g-1 
                   
                   if (reference@data@withContinuous){
                     rownames(reference@param@paramContinuous@mu)  <-   colnames(reference@data@dataContinuous@data)
                     rownames(reference@param@paramContinuous@sd)  <-   colnames(reference@data@dataContinuous@data)
                     reference@param@paramContinuous@pi <-  numeric()
                     colnames(reference@param@paramContinuous@mu) <-  paste("class-",1:length(reference@param@pi),sep="")
                     colnames(reference@param@paramContinuous@sd) <-  paste("class-",1:length(reference@param@pi),sep="")
                     shortomega <- reference@model@omega[which(names(reference@model@omega) %in% colnames(reference@data@dataContinuous@data))]
                     nbparam <- nbparam + sum(shortomega) * reference@model@g * 2 + sum(1-shortomega) * 2
                   }
                   
                   if (reference@data@withInteger){
                     rownames(reference@param@paramInteger@lambda)  <-   colnames(reference@data@dataInteger@data)
                     reference@param@paramInteger@pi <- numeric()
                     colnames(reference@param@paramInteger@lambda) <-  paste("class-",1:length(reference@param@pi),sep="")
                     shortomega <- reference@model@omega[which(names(reference@model@omega) %in% colnames(reference@data@dataInteger@data))]
                     nbparam <- nbparam + sum(shortomega) * reference@model@g  + sum(1-shortomega)
                   }
                   if (reference@data@withCategorical){
                     reference@param@paramCategorical@pi <- numeric()
                     for (j in 1:reference@data@dataCategorical@d){
                       reference@param@paramCategorical@alpha[[j]] <- matrix(reference@param@paramCategorical@alpha[[j]], nrow = reference@model@g)
                       rownames(reference@param@paramCategorical@alpha[[j]]) <- paste("class-",1:length(reference@param@pi),sep="")
                       colnames(reference@param@paramCategorical@alpha[[j]]) <- reference@data@dataCategorical@modalitynames[[j]]
                     }                     
                     names(reference@param@paramCategorical@alpha) <- colnames(reference@data@dataCategorical@shortdata)
                     shortomega <- reference@model@omega[which(names(reference@model@omega) %in% colnames(reference@data@dataCategorical@shortdata))]
                     for (j in 1:length(shortomega)){
                       nbparam <- nbparam + (length(reference@data@dataCategorical@modalitynames[[j]])-1) * (1 + (reference@model@g-1)*shortomega[j])
                     }
                   }
                   reference@criteria@nbparam <- nbparam
                   reference@criteria@BIC <- reference@criteria@loglikelihood - nbparam*0.5*log(reference@data@n)
                   reference@criteria@AIC <- reference@criteria@loglikelihood - nbparam
                   names(reference@criteria@BIC) <- NULL
                   reference@criteria@ICL <- ICLmixed(reference) 
                 }else{
                   warning("All the models got error (degeneracy)", call. = FALSE)  
                 }
               }
             }
             if (reference@data@withContinuous){
               for (j in 1:reference@data@dataContinuous@d){
                 if (any(reference@data@dataContinuous@notNA[,j]==0))
                   reference@data@dataContinuous@data[which(reference@data@dataContinuous@notNA[,j]==0),j] <- NA
               }
             }
             if (reference@data@withInteger){
               for (j in 1:reference@data@dataInteger@d){
                 if (any(reference@data@dataInteger@notNA[,j]==0))
                   reference@data@dataInteger@data[which(reference@data@dataInteger@notNA[,j]==0),j] <- NA
               }
             }
             if (reference@data@withCategorical){
               if (any(reference@data@dataCategorical@shortdata == 0)){
                 for (j in 1:ncol(reference@data@dataCategorical@shortdata))
                   reference@data@dataCategorical@shortdata[which(reference@data@dataCategorical@shortdata[,j] == 0), ] <- NA
               }
             }
             
             names(reference@criteria@nbparam) <- NULL
             if (any(reference@model@omega==1))
               reference@model@names.relevant <- as.character(names(reference@model@omega)[which(reference@model@omega==1)])
             if (any(reference@model@omega==0))
               reference@model@names.irrelevant <- names(reference@model@omega)[which(reference@model@omega==0)]
             if (reference@strategy@crit.varsel != "MICL") reference@criteria@MICL <- numeric()
             return(reference)
           }
)