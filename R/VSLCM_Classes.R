setClass(
  Class = "VSLCMcriteria", 
  representation = representation(
    likelihood="numeric",
    BIC="numeric",
    ICL="numeric",
    MICL="numeric"
  ), 
  prototype = prototype(
    likelihood=numeric(),
    BIC=numeric(),
    ICL=numeric(),
    MICL=numeric()
  )
)

setClass(
  Class = "VSLCMpartitions", 
  representation = representation(
    zMAP="numeric",
    zOPT="numeric",
    tik="matrix"
  ), 
  prototype = prototype(
    zMAP=numeric(),
    zOPT=numeric(),
    tik=matrix(0,0,0)
  )
)




setClass(
  Class = "VSLCMmodel", 
  representation = representation(
    g="numeric",
    omega="numeric"
  ), 
  prototype = prototype(
    g=numeric(),
    omega=numeric()
  )
)


setClass(
  Class = "VSLCMparameters", 
  representation = representation(
    proportions="numeric",
    means="matrix",
    variances="matrix"
  ), 
  prototype = prototype(
    proportions=numeric(),
    means=matrix(),
    variances=matrix()
  )
)



setClass(
  Class = "VSLCMresults", 
  representation = representation(
    data="matrix",
    priors="matrix",
    criteria="VSLCMcriteria",
    partitions="VSLCMpartitions",
    model="VSLCMmodel",
    parameters="VSLCMparameters"
  ), 
  prototype = prototype(
    data=matrix(),
    priors=matrix(),
    criteria=new("VSLCMcriteria"),
    partitions=new("VSLCMpartitions"),
    model=new("VSLCMmodel"),
    parameters=new("VSLCMparameters")
  )
)







