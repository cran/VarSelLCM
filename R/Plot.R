

plotCateg <- function(data, model, param){
  
  results <- matrix(0,model@g,data@d)
  for (j in 1:data@d){
    for (k in 1:model@g){
      alpha <- 0
      for (k2 in (c(1:model@g)[-k]))
        alpha <- alpha + (param@pi[k2]* param@alpha[[j]][k2,])/sum(param@pi[-k])
      
      results[k,j] <- round(sum((param@alpha[[j]][k,]-alpha)**2),4)
    }
  }
  
  for (j in  which(model@omega==1)){
  mp<-barplot(matrix(results[,j], ncol = 1), beside = TRUE, axisnames = FALSE, main=(colnames(data@shortdata))[j],  ylab="Discriminative measure", ylim=c(0,max(c(results,1))), cex.main=0.8)
  mtext(1, at = mp, text = paste("C",1:model@g), line = 0, cex = 0.5)
  }
}

setMethod(
  f="plot",
  signature = c("VSLCMresultsCategorical"),
  definition = function(x){
    op <- par(no.readonly = TRUE)
    if (any(x@model@omega==1) && (x@model@g>1)){
      nvar <- sum(x@model@omega)
      # split the layout
      if( nvar < 4 & nvar > 1){
        par( mfrow = c( 1, nvar ), mar=c(2,4,2,2), cex.axis = 0.6, cex.lab = 0.6, cex.main = 0.7  )
      }else if ( nvar >= 4 ){
        nrow<-round(sqrt(nvar))
        ncol <- ceiling(nvar/nrow)
        par( mfrow = c( nrow, ncol ), mar=c(2,4,2,2), cex.axis = 0.6, cex.lab = 0.6, cex.main = 0.7  ) 
      }
      plotCateg(x@data, x@model, x@param)
    }
    else
      cat("No plot is available since none variable is discriminative!")
    par(op)
  }
  
)

plotCont <- function(data, model, param){
  
  results <- matrix(0,model@g,data@d)
  for (j in 1:data@d){
    canddef <- range(data@data[,j], na.rm = T)
    canddef <- canddef + c(-0.05,0.05)*(canddef[2] - canddef[1])
    for (k in 1:model@g){
      ftmp <- function(u){
        out <- dnorm(u, param@mu[j,k], param@sd[j,k])
        for (k2 in(c(1:model@g)[-k])) out <- out - param@pi[k2]/sum(param@pi[-k]) * dnorm(u, param@mu[j,k2], param@sd[j,k2])
        return(out**2)
      }
      results[k,j] <- round(integrate(ftmp, lower = canddef[1], upper = canddef[2], subdivisions = 10**5)$value,4)
    }
  }
  
  for (j in which(model@omega==1)){
    mp<-barplot(matrix(results[,j], ncol = 1), beside = TRUE, axisnames = FALSE, main=(colnames(data@data))[j], ylab="Discriminative measure", ylim=c(0,max(c(results,1))), cex.main=0.8)
    mtext(1, at = mp, text = paste("C",1:model@g), line = 0, cex = 0.5)
  }
}

setMethod(
  f="plot",
  signature = c("VSLCMresultsContinuous"),
  definition = function(x){
    op <- par(no.readonly = TRUE)
    if (any(x@model@omega==1) && (x@model@g>1)){
      nvar <- sum(x@model@omega)
      # split the layout
      if( nvar < 4 & nvar > 1){
        par( mfrow = c( 1, nvar ), mar=c(2,4,2,2), cex.axis = 0.6, cex.lab = 0.6, cex.main = 0.7  )
      }else if ( nvar >= 4 ){
        nrow<-round(sqrt(nvar))
        ncol <- ceiling(nvar/nrow)
        par( mfrow = c( nrow, ncol ), mar=c(2,4,2,2), cex.axis = 0.6, cex.lab = 0.6, cex.main = 0.7  ) 
      }
      plotCont(x@data, x@model, x@param)    
    }
    else
      cat("No plot is available since none variable is discriminative!")
    par(op)
  }
)


plotInte <- function(data, model, param){
  results <- matrix(0,model@g,data@d)
  for (j in 1:data@d){
    canddef <- range(data@data[,j], na.rm = T)
    for (k in 1:model@g){
      ftmp <- function(u){
        out <- dpois(u, param@lambda[j,k])
        for (k2 in(c(1:model@g)[-k])) out <- out - param@pi[k2]/sum(param@pi[-k]) * dpois(u, param@lambda[j,k2])
        return(out**2)
      }
      results[k,j] <- round(mean(ftmp(canddef[1] : canddef[2])),4)
    }
  }
  for (j in which(model@omega==1)){
    mp<-barplot(matrix(results[,j], ncol = 1), beside = TRUE, axisnames = FALSE,  main=(colnames(data@data))[j], ylab="Discriminative measure", ylim=c(0,max(c(results,1))), cex.main=0.8)
    mtext(1, at = mp, text = paste("C",1:model@g), line = 0, cex = 0.5)
  }
  
}

setMethod(
  f="plot",
  signature = c("VSLCMresultsInteger"),
  definition = function(x){
    op <- par(no.readonly = TRUE)
    if (any(x@model@omega==1) && (x@model@g>1)){
      nvar <- sum(x@model@omega)
      # split the layout
      if( nvar < 4 & nvar > 1){
        par( mfrow = c( 1, nvar ), mar=c(2,4,2,2), cex.axis = 0.6, cex.lab = 0.6, cex.main = 0.7  )
      }else if ( nvar >= 4 ){
        nrow<-round(sqrt(nvar))
        ncol <- ceiling(nvar/nrow)
        par( mfrow = c( nrow, ncol ), mar=c(2,4,2,2), cex.axis = 0.6, cex.lab = 0.6, cex.main = 0.7  ) 
      }
      plotInte(x@data, x@model, x@param)    
    }
    else
      cat("No plot is available since none variable is discriminative!")
    par(op)
  }
)


setMethod(
  f="plot",
  signature = c("VSLCMresultsMixed"),
  definition = function(x){
    op <- par(no.readonly = TRUE)
    if ((x@model@g>1) && (sum(x@model@omega)>0)){
      nvar <- sum(x@model@omega)
      # split the layout
      if( nvar < 4 & nvar > 1){
        par( mfrow = c( 1, nvar ), mar=c(2,4,2,2), cex.axis = 0.6, cex.lab = 0.6, cex.main = 0.7  )
      }else if ( nvar >= 4 ){
        nrow<-round(sqrt(nvar))
        ncol <- ceiling(nvar/nrow)
        par( mfrow = c( nrow, ncol ), mar=c(2,4,2,2), cex.axis = 0.6, cex.lab = 0.6, cex.main = 0.7 ) 
      }
      
      if (x@data@withContinuous)
        plotCont(x@data@dataContinuous, new("VSLCMmodel", g=x@model@g, omega=x@model@omega[which(names(x@model@omega) %in% colnames(x@data@dataContinuous@data))]), x@param@paramContinuous)    
      
      if (x@data@withContinuous)
        plotInte(x@data@dataInteger, new("VSLCMmodel", g=x@model@g, omega=x@model@omega[which(names(x@model@omega) %in% colnames(x@data@dataInteger@data))]), x@param@paramInteger)    
      
      if (x@data@withCategorical)
        plotCateg(x@data@dataCategorical, new("VSLCMmodel", g=x@model@g, omega=x@model@omega[which(names(x@model@omega) %in% colnames(x@data@dataCategorical@shortdata))]), x@param@paramCategorical)    
      
    }
    par(op)
  }
)

