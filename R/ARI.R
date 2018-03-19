###################################################################################
##' Adjusted Rand Index
##'
##' @description  
##' This function computes the Adjusted Rand Index
##'
##' @param x vector defining a partition.
##' @param y vector defining a partition of whose length is equal to the length of x.
##' 
##' @return numeric
##' 
##' @references L. Hubert and P. Arabie (1985) Comparing Partitions, Journal of the Classification, 2, pp. 193-218. 
##' @examples
##' x <- sample(1:2, 20, replace=TRUE)
##' y <- x
##' y[1:5] <- sample(1:2, 5, replace=TRUE)
##' ARI(x, y)
##' @export
##'
##'
ARI <- function (x, y) 
{
  if ((length(x) != length(y)))
    stop("The two partitions must be vectors of same length")
  ari <- 0
  if ((length(unique(x)) + length(unique(y))) == 2) ari <- 1
  conting <- table(x, y)
  a <- sum(choose(conting, 2))
  b <- sum(choose(rowSums(conting), 2)) - a
  c <- sum(choose(colSums(conting), 2)) - a
  d <- choose(sum(conting), 2) - a - b - c
  ari <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + 
                                                     a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ari)
}