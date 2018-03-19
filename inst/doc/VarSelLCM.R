## ---- comment=""---------------------------------------------------------
library(VarSelLCM)
# Data loading
data("heart")
head(heart)

## ---- comment=""---------------------------------------------------------
# Add a missing value artificially (just to show that it works!)
heart[1,1] <- NA
# Clustering with variable selection and a number of cluster betwee 1 and 4
# Model selection is BIC (to use MICL, the option must be specified)
out <- VarSelCluster(heart[,-13], 1:4, nbcores = 2)

## ---- eval=FALSE, comment="", include=TRUE-------------------------------
#  # Start the shiny application
#   VarSelShiny(out)

## ---- comment=""---------------------------------------------------------
# Summary of the best model
summary(out)

## ---- comment=""---------------------------------------------------------
# Discriminative power of the variables (here, the most discriminative variable is MaxHeartRate)
plot(out, type="bar")

## ---- comment=""---------------------------------------------------------
# Boxplot for continuous (or interger) variable
plot(out, y="MaxHeartRate", type="boxplot")

## ---- comment=""---------------------------------------------------------
# Empirical and theoretical distributions (to check that clustering is pertinent)
plot(out, y="MaxHeartRate", type="cdf")

## ---- comment=""---------------------------------------------------------
# Summary of categorical variable
plot(out, y="Sex")

## ---- comment=""---------------------------------------------------------
# Summary of the probabilities of missclassification
plot(out, type="probs-class")

## ---- comment=""---------------------------------------------------------
# Imputation by posterior mean for the first observation
not.imputed <- heart[1,-13]
imputed <- VarSelImputation(out)[1,]
rbind(not.imputed, imputed)

