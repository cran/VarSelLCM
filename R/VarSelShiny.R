###################################################################################
###################################################################################
##' Shiny app for analyzing results from VarSelCluster
##'
##' @description  
##' Shiny app for analyzing results from VarSelCluster
##' 
##' @param X an instance of \linkS4class{VSLCMresults} returned by function \link{VarSelCluster}.
##' 
##' @examples
##' \dontrun{
##' # Data loading
##' data("heart")
##' # Clustering en 2 classes
##' results <- VarSelCluster(heart[,-13], 2)
##' # Opening Shiny application to easily see the results
##' VarSelShiny(results)
##' }
##' 
##' @export
##'
##'
VarSelShiny <-  function(X){
  check.results(X)
  G <- .GlobalEnv
  assign("resVSLC", X, envir=G)
  a=shiny::runApp(system.file("shinyApp",package="VarSelLCM"),launch.browser = TRUE)
  return(invisible(a))
}
