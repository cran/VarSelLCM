# Verifie les parametres d entrees
CheckInputs <- function(x, g, initModel, vbleSelec, crit.varsel, discrim, paramEstim, nbcores, nbSmall, iterSmall, nbKeep, iterKeep, tolKeep){
  if ( (is.numeric(g)==FALSE) || (length(g)!=1))
    stop("The component number have to be an integer of length one!")
  
  if (is.matrix(x))
    x <- as.data.frame(x)
  
  if (is.data.frame(x)==FALSE)
    stop("Data set must be a data frame!")
  
  if (is.logical(vbleSelec) == FALSE)
    stop("Input vbleSelec must be logicial")
  
  if (!(crit.varsel %in% c("AIC", "BIC", "ICL")) & (!vbleSelec))
    stop("Input vbleSelec must be equal to AIC, BIC or ICL without variable selection")
  
  if (!(crit.varsel %in% c("AIC", "BIC", "MICL")) & (vbleSelec))
    stop("Input vbleSelec must be equal to AIC, BIC or MICL with variable selection")
  
  if ((length(discrim) != ncol(x)) || (all(discrim %in% c(0,1))==FALSE))
    stop("Input discrim must be logical of length number of variables")
  
  if (is.logical(paramEstim) == FALSE)
    stop("Input paramEstim must be logicial")
  
  if (is.numeric(nbcores) == FALSE)
    stop("Input nbcores must be numeric")
  
  if ((is.numeric(nbSmall) == FALSE) || (length(nbSmall)!=1))
    stop("Input nbSmall must be numeric of size one")
  
  if ((is.numeric(iterSmall) == FALSE) || (length(iterSmall)!=1))
    stop("Input iterSmall must be numeric of size one")
  
  if ((is.numeric(nbKeep) == FALSE) || (length(nbKeep)!=1))
    stop("Input nbKeep must be numeric of size one")
  
  if ((is.numeric(iterKeep) == FALSE) || (length(iterKeep)!=1))
    stop("Input iterKeep must be numeric of size one")
  
  if ((is.numeric(tolKeep) == FALSE) || (length(tolKeep)!=1))
    stop("Input tolKeep must be numeric of size one")  
}