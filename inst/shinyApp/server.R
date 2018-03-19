
shinyServer(function(input, output) {
  output$summary <- renderPrint(summary(resVSLC))
  
  output$plot1 <- renderPlot(plot(resVSLC, type=input$typeDiscrim, ylim=c(1, input$numDiscrim)))

  output$table1 <- renderTable({
    info <- resVSLC@criteria@discrim[which(resVSLC@criteria@discrim>0)]
    tmp <- data.frame(
      variables=names(info),
      discrim=info,
      discrim.percent=100*info/sum(info) ,
      discrim.cumsum=100*cumsum(info)/sum(info)) 
    colnames(tmp) <- c("Variables", "Discrim. Power", "Discrim. Power (%)", "Discrim. Power (% cum)")
    tmp[1:input$numDiscrim,]
  }
  )

  output$plot2 <- renderPlot( plot(resVSLC, y=input$vble, type=input$typevble)) 
  
  output$plot3 <- renderPlot(plot(x=resVSLC, type=input$probs))
})
