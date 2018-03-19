########################################################################################################################
## Fonction calculant la vraisemblance completee integree en zMAP (ICL exact)
########################################################################################################################
## fonction outils
logdinvgamma <- function(x, alpha, beta)
  alpha*log(beta) - lgamma(alpha) - (alpha+1)*log(x) - beta/x

## Integrale sur une variable continue pour 1 classe ou non discriminante
IntegreOneVariableContinuous <- function(x, priors){
  nu <- priors[1]
  s0 <- priors[2]
  mu <- priors[3]
  n0 <- priors[4]
  n <- length(x)
  n1 <- n + n0
  s1 <- sqrt( sum((x-mean(x))**2) + s0**2 + ((mu - mean(x))**2)/(1/n0 + 1/n) )
  -log(sqrt(pi))*n + lgamma((n + nu)/2) - lgamma(nu/2) +   nu * log(s0/s1) - n*log(s1) + log(sqrt(n0 / n1) )

}
## Integrale sur une variable integer pour 1 classe ou non discriminante
IntegreOneVariableInteger <- function(x, priors){
  alpha <- sum(x) + priors[1]
  beta <- length(x) + priors[2]
  integre <- priors[1] * log(priors[2]) - lgamma(priors[1]) + lgamma(alpha) - alpha * log(beta) - sum(lgamma(x+1))
  return(integre)
}

## Integrale sur une variable categorielle pour 1 classe ou non discriminante
IntegreOneVariableCategorical <- function(eff) lgamma(length(eff) * 0.5) - length(eff) * lgamma(0.5) + sum(lgamma(eff + 0.5)) - lgamma(sum(eff + 0.5))

## ICL exact dans le cas de variables mixed
pvdiscrim <- function(obj){
  out <- rep(0, obj@data@d)
  repere <- 0
  if (obj@data@withContinuous){
    for (j in 1:obj@data@dataContinuous@d){
      repere <- repere + 1
      out[repere] <- sum(sapply(1:obj@model@g, function(k) IntegreOneVariableContinuous(obj@data@dataContinuous@data[which( (obj@partitions@zMAP==k)*obj@data@dataContinuous@notNA[,j] ==1),j], obj@data@dataContinuous@priors[j,]))) - IntegreOneVariableContinuous(obj@data@dataContinuous@data[which(obj@data@dataContinuous@notNA[,j]==1),j], obj@data@dataContinuous@priors[j,])
      names(out)[repere] <- colnames(obj@data@dataContinuous@data)[j]
    }
  }
  if (obj@data@withInteger){
    for (j in 1:obj@data@dataInteger@d){
      repere <- repere + 1
      out[repere] <- sum(sapply(1:obj@model@g, function(k) IntegreOneVariableInteger(obj@data@dataInteger@data[which( (obj@partitions@zMAP==k)*obj@data@dataInteger@notNA[,j] ==1),j], obj@data@dataInteger@priors[j,]))) - IntegreOneVariableInteger(obj@data@dataInteger@data[which(obj@data@dataInteger@notNA[,j]==1),j], obj@data@dataInteger@priors[j,])
      names(out)[repere] <- colnames(obj@data@dataInteger@data)[j]
    }
  }
  if (obj@data@withCategorical){
    who <- which(names(obj@model@omega) %in% colnames(obj@data@dataCategorical@shortdata))
    for (loc in 1:obj@data@dataCategorical@d){
      repere <- repere + 1
      names(out)[repere] <- colnames(obj@data@dataCategorical@data)[loc]
      eff <- table(obj@data@dataCategorical@data[,loc])
      out[repere] <-  -IntegreOneVariableCategorical(eff)
      
      for (k in 1:obj@model@g){
        for (h in 1:length(eff))  eff[h] <- sum(obj@data@dataCategorical@data[which(obj@partitions@zMAP == k),loc] == obj@data@dataCategorical@modalitynames[[loc]][h], na.rm = T)
        out[repere] <- out[repere]  + IntegreOneVariableCategorical(eff)
      }
    }    
  }
  return(out)
}



## ICL exact dans le cas de variables continues
ICLcontinuous <- function(obj){
  ICLexact <- lgamma(obj@model@g/2) - obj@model@g*lgamma(1/2) + 
    sum(lgamma(table(c(1:obj@model@g, obj@partitions@zMAP)) - 1/2)) - lgamma(length(obj@partitions@zMAP) + obj@model@g/2)
  for (j in 1:obj@data@d){
    if (obj@model@omega[j]==0){
      ICLexact <- ICLexact  + IntegreOneVariableContinuous(obj@data@data[which(obj@data@notNA[,j]==1),j], obj@data@priors[j,])
    }else{
      for (k in unique(obj@partitions@zMAP))
        ICLexact <- ICLexact  + IntegreOneVariableContinuous(obj@data@data[which( (obj@partitions@zMAP==k)*obj@data@notNA[,j] ==1),j], obj@data@priors[j,])
    }
  }
  names(ICLexact) <- NULL
  return(ICLexact)
}

## ICL exact dans le cas de variables entiere
ICLinteger <- function(obj){
  ICLexact <- lgamma(obj@model@g/2) - obj@model@g*lgamma(1/2) + 
    sum(lgamma(table(c(1:obj@model@g, obj@partitions@zMAP)) - 1/2)) - lgamma(length(obj@partitions@zMAP) + obj@model@g/2)
  for (j in 1:obj@data@d){
    if (obj@model@omega[j]==0){
      ICLexact <- ICLexact  + IntegreOneVariableInteger(obj@data@data[which(obj@data@notNA[,j]==1),j], obj@data@priors[j,])
    }else{
      for (k in unique(obj@partitions@zMAP))
        ICLexact <- ICLexact  + IntegreOneVariableInteger(obj@data@data[which( (obj@partitions@zMAP==k)*obj@data@notNA[,j] ==1),j], obj@data@priors[j,])
    }
  }
  names(ICLexact) <- NULL
  return(ICLexact)
}


## ICL exact dans le cas de variables categorielles
ICLcategorical <- function(obj){
  ICLexact <- lgamma(obj@model@g/2) - obj@model@g*lgamma(1/2) + sum(lgamma(table(c(1:obj@model@g, obj@partitions@zMAP)) - 1/2)) - lgamma(obj@data@n + obj@model@g/2)
  for (j in 1:obj@data@d){
    eff <- rep(0, length(obj@data@modalitynames[[j]]))
    if (obj@model@omega[j]==0){
      for (h in 1:length(eff))
        eff[h] <- sum(obj@data@data[,j] == obj@data@modalitynames[[j]][h], na.rm = TRUE)
      ICLexact <- ICLexact  + IntegreOneVariableCategorical(eff)
    }else{
      for (k in 1:obj@model@g){
        for (h in 1:length(eff))
          eff[h] <- sum(obj@data@data[which(obj@partitions@zMAP == k),j] == obj@data@modalitynames[[j]][h], na.rm = TRUE)
        ICLexact <- ICLexact  + IntegreOneVariableCategorical(eff)        
      }
    }
  }
  names(ICLexact) <- NULL
  return(ICLexact)
}


## ICL exact dans le cas de variables mixed
ICLmixed <- function(obj){
  ICLexact <- lgamma(obj@model@g/2) - obj@model@g*lgamma(1/2) +  sum(lgamma(table(c(1:obj@model@g, obj@partitions@zMAP)) - 1/2)) - lgamma(length(obj@partitions@zMAP) + obj@model@g/2)
  if (obj@data@withContinuous){
    who <- which(names(obj@model@omega) %in% colnames(obj@data@dataContinuous@data))
    loc <- 0
    for (j in who){
      loc <- loc + 1
      if (obj@model@omega[j]==0){
        ICLexact <- ICLexact  + IntegreOneVariableContinuous(obj@data@dataContinuous@data[which(obj@data@dataContinuous@notNA[,loc]==1),loc], obj@data@dataContinuous@priors[loc,])
      }else{
        for (k in unique(obj@partitions@zMAP))
          ICLexact <- ICLexact  + IntegreOneVariableContinuous(obj@data@dataContinuous@data[which( (obj@partitions@zMAP==k)*obj@data@dataContinuous@notNA[,loc] ==1),loc], obj@data@dataContinuous@priors[loc,])
      }
    }
  }
  
  if (obj@data@withInteger){
    who <- which(names(obj@model@omega) %in% colnames(obj@data@dataInteger@data))
    loc <- 0
    for (j in who){
      loc <- loc + 1
      if (obj@model@omega[j]==0){
        ICLexact <- ICLexact  + IntegreOneVariableInteger(obj@data@dataInteger@data[which(obj@data@dataInteger@notNA[,loc]==1),loc], obj@data@dataInteger@priors[loc,])
      }else{
        for (k in unique(obj@partitions@zMAP))
          ICLexact <- ICLexact  + IntegreOneVariableInteger(obj@data@dataInteger@data[which( (obj@partitions@zMAP==k)*obj@data@dataInteger@notNA[,loc] ==1),loc], obj@data@dataInteger@priors[loc,])
      }
    }
  }
  
  if (obj@data@withCategorical){
    who <- which(names(obj@model@omega) %in% colnames(obj@data@dataCategorical@shortdata))
    loc <- 0
    for (j in who){
      loc <- loc + 1
      eff <- rep(0, length(obj@data@dataCategorical@modalitynames[[loc]]))
      if (obj@model@omega[j]==0){
        for (h in 1:length(eff))
          eff[h] <- sum(obj@data@dataCategorical@data[,loc] == obj@data@dataCategorical@modalitynames[[loc]][h], na.rm = TRUE)
        ICLexact <- ICLexact  + IntegreOneVariableCategorical(eff)
      }else{
        for (k in unique(obj@partitions@zMAP)){
          for (h in 1:length(eff))
            eff[h] <- sum(obj@data@dataCategorical@data[which(obj@partitions@zMAP == k),loc] == obj@data@dataCategorical@modalitynames[[loc]][h], na.rm = T)
          ICLexact <- ICLexact  + IntegreOneVariableCategorical(eff)
        }
      }    
    }
  }
  names(ICLexact) <- NULL
  return(ICLexact)
}

