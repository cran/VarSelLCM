shinyUI(fluidPage(
  titlePanel("Results of the R package VarSelLCM"),
  
  sidebarPanel(
    h3("Theoretical details"),
    tags$div(class="header", checked=NA, tags$a(href="https://link.springer.com/article/10.1007/s11222-016-9670-1", "Read the paper")),
   h3("Discrimintative power"),
   radioButtons("typeDiscrim", "Graph", c("Bar-chart" = "bar", "Pie" = "pie")),
   sliderInput("numDiscrim", "Number of the most discriminative variables:", min =1, max = sum(resVSLC@model@omega), value = sum(resVSLC@model@omega)),
   h3("Univariate Analysis"),
   selectizeInput("vble", "Name of the variable", resVSLC@data@var.names, selected = NULL, multiple = FALSE),
   radioButtons("typevble", "Type", c("CDF" = "cdf", "Boxplot" = "boxplot")),
   h3("Misclassification probabilities"),
   radioButtons("probs", "Histrogram", c("Overall" = "probs-overall", "Per class" = "probs-class"))
  ),
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Summary", verbatimTextOutput("summary")),
                tabPanel("Discriminative power", plotOutput("plot1"), tableOutput("table1")),
                tabPanel("Univariate Analysis", plotOutput("plot2")),
                tabPanel("Misclassification probabilities", plotOutput("plot3"))
    )   
  )
))
