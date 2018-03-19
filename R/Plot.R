cdfmixtureGauss <- function(u, pi, mu, sd){
  out <- u * 0
  for (k in 1:length(pi)) 
    out <- out + pnorm(u, mu[k], sd[k])*pi[k]
  out
}

varsellcm.plot.cont.cdf <- function(df, y, param){
  graphic <- ggplot(df, aes(x=df$x)) +   
    stat_ecdf(geom = "step", aes(colour="Empirical"), size=1) +
    stat_function(fun = cdfmixtureGauss, args = param, aes(colour = "Theoretical"), size=1) +
    scale_colour_manual("CDF", values = c("black", "dodgerblue3")) +
    scale_x_continuous(name = y) + 
    scale_y_continuous(name="CDF") + 
    ggtitle(paste("CDF of", y)) +  
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  print(graphic)
}

varsellcm.plot.boxplot <- function(df, y){
  graphic <- ggplot(df, aes(x=class, y=df$x, fill=class)) + 
    geom_boxplot() +   
    guides(fill=FALSE) +
    coord_flip() +
    scale_y_continuous(name=y)  + 
    ggtitle(paste("Boxplots of", y)) +  
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  print(graphic)
}

cdfmixturePoiss <- function(u, pi, lam){
  out <- u * 0
  for (k in 1:length(pi)) 
    out <- out + ppois(u, lam[k])*pi[k]
  out
}

varsellcm.plot.inte.cdf <- function(df, y, param){
  graphic <- ggplot(df, aes(x=df$x)) +   
    stat_ecdf(geom = "step", aes(colour="Empirical"), size=1) +
    stat_function(fun = cdfmixturePoiss, args = param, aes(colour = "Theoretical"), size=1) +
    scale_colour_manual("CDF", values = c("black", "dodgerblue3")) +
    scale_x_continuous(name = y) + 
    scale_y_continuous(name="CDF") + 
    ggtitle(paste("CDF of", y)) +  
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  print(graphic)
}

varsellcm.plot.cate  <- function(tmp, y){
  
  df <- data.frame(class=as.factor(rep(1:nrow(tmp), ncol(tmp))),
                   levels=rep(colnames(tmp), each=nrow(tmp)),
                   probabilties=round(as.numeric(tmp), 2))
  
  
  
  graph <- ggplot(data=df, aes(x=levels, y=df$probabilties, fill=class)) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(aes(label=df$probabilties), vjust=1.6, color="black",  position = position_dodge(0.9), size=3.5)+
    scale_fill_brewer(palette="Paired")+
    scale_y_continuous(name="probability")+
    theme_minimal() +
    ggtitle(paste("Distribution per class of", y)) +  
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  print(graph)
}


#' This function draws information about an instance of \code{\linkS4class{VSLCMresults}}.
#' 
#' @param x instance of  \code{\linkS4class{VSLCMresults}}.
#' @param y character. The name of the variable we are interested in. Otherwise, graphic about general results is ploted. 
#' @param type character. Can be cdf or boxplot to plot the distribution of the variable we are interested in. Otherwise, must be pie or bar to obtain information about the discriminative power of the variables and must be probs-overall or probs-class to obtain information about the misclassification probabilities. 
#' @param ylim numeric. Define the range of the most discriminative variables to considered (only use if type="pie" or type="bar")
#' @param ... Additional argument list that might not ever be used.
#' 
#' @name plot
#' @rdname plot-methods
#' @docType methods
#' @exportMethod plot
NULL
#' @rdname plot-methods
#' @aliases plot plot,VSLCMresults,ANY-method
#' @aliases plot plot,VSLCMresults,character,ANY-method
#' 
setMethod(f="plot",
  signature = c("VSLCMresults", "character"),
  definition = function(x, y, type,ylim=NULL){
    vu <- FALSE
    if (x@data@withContinuous){
      if (y %in% rownames(x@param@paramContinuous@mu)){
        loc2 <- which(rownames(x@param@paramContinuous@mu)==y)
        if (length(loc2)!=1)
          stop("y must be the name of a variable in the analyzed data")
        if (type=="cdf")
          varsellcm.plot.cont.cdf(data.frame(x = x@data@dataContinuous@data[, which(rownames(x@param@paramContinuous@mu)==y)]),
                                  y,
                                  list(x@param@pi, x@param@paramContinuous@mu[loc2,], x@param@paramContinuous@sd[loc2,]))
        else if (type=="boxplot")
          varsellcm.plot.boxplot(data.frame(x = x@data@dataContinuous@data[, which(rownames(x@param@paramContinuous@mu)==y)],
                                            class=as.factor(x@partitions@zMAP)),
                                 y)
        else
          stop("type must be cdf or boxplot")
        vu <- TRUE
      }
    }
    if (x@data@withInteger){
      if (y %in% rownames(x@param@paramInteger@lambda)){
        loc2 <- which(rownames(x@param@paramInteger@lambda)==y)
        if (length(loc2)!=1)
          stop("y must be the name of a variable in the analyzed data")
        if (type=="cdf")
          varsellcm.plot.inte.cdf(data.frame(x = x@data@dataInteger@data[, which(rownames(x@param@paramInteger@lambda)==y)]),
                                  y,
                                  list(x@param@pi, x@param@paramInteger@lambda[loc2,]))
        else if (type=="boxplot")
          varsellcm.plot.boxplot(data.frame(x = x@data@dataInteger@data[, which(rownames(x@param@paramInteger@lambda)==y)],
                                            class=as.factor(x@partitions@zMAP)),
                                 y)
        else
          stop("type must be cdf or boxplot")
        vu <- TRUE
      }
    }
    if (x@data@withCategorical){
      if (y %in% names(x@param@paramCategorical@alpha)){
        loc2 <- which(names(x@param@paramCategorical@alpha) ==y)
        ifelse (length(loc2)==1,
                varsellcm.plot.cate(x@param@paramCategorical@alpha[[loc2]], y),
                stop("y must be the name of a variable in the analyzed data"))

        vu <- TRUE
      }
    }
    if (!vu)
      stop("y must be the name of a variable in the analyzed data")
  }
)



setMethod(
  f="plot",
  signature = c("VSLCMresults"),
  definition = function(x, type, ylim=c(1, x@data@d)){
    df <- data.frame(discrim.power=x@criteria@discrim, variables=as.factor(names(x@criteria@discrim)), rg=1:x@data@d)
    df <- df[which(df$discrim.power>0),]
    df <- df[order(df$discrim.power, decreasing = T),]
    ylim <- as.integer(sort(ylim))
    ylim[2] <- min(ylim[2], nrow(df))
    ylim[1] <- min(max(ylim[1], 0), ylim[2])
    df <- df[ylim[1]:ylim[2], , drop=F]
    if (type=="pie"){
      pie<- ggplot(df, aes(x="", y=df$discrim.power, fill=df$variables))+
        scale_y_continuous(name="discriminative power") +
        geom_bar(width = 1, stat = "identity") +
        coord_polar("y", start=0)  +
        ggtitle(paste("Discriminative power")) +
        theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
      print(pie)
    }else if (type=="bar"){
      bar <- ggplot(data=df, aes(x=df$rg, y=df$discrim.power, fill=df$variables)) +
        scale_y_continuous(name="discriminative power") +
        geom_bar(stat="identity", position=position_dodge())+
        geom_text(aes(label=round(df$discrim.power,2)), vjust=-0.1, color="black",
                  position = position_dodge(0.9), size=3.5)+
        scale_x_discrete(name="Variables")+
        scale_fill_brewer(palette="Paired")+
        theme_minimal()  +
        ggtitle(paste("Discriminative power")) +
        theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
      print(bar)
    }else if (type=="probs-overall"){
      tmp <- data.frame(probs=1-apply(x@partitions@tik, 1, max))
      tikplot <- ggplot(tmp, aes(tmp$probs)) +   geom_histogram(binwidth = 0.05) + scale_x_continuous("Probability of misclassification") +
        ggtitle(paste("Probabilities of misclassification")) +
        theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
      print(tikplot)
    }else if (type=="probs-class"){
      tmp <- data.frame(probs=1-apply(x@partitions@tik, 1, max), class=as.factor(x@partitions@zMAP))
      tikplot <-    ggplot(tmp, aes(x=tmp$probs, fill=class)) +   geom_histogram(position="dodge", binwidth = 0.05)+ scale_x_continuous("Probability of misclassification") +
        ggtitle(paste("Probabilities of misclassification")) +
        theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
      print(tikplot)
    }else{
      stop("type must be specified and equal to pie or bar or probs-overall or class")
    }
  }
)