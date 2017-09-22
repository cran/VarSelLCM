\name{heart}
\alias{heart}
\docType{data}
\title{
Heart desease data set
}
\description{
Heart data set composed with 270 individuals described by three continuous variables, one count variable and eight categorical variables. Individuals are split into two groups.
}
\usage{data("heart")}
\format{
  A data frame with 270 observations on the following 13 variables.
  \describe{
    \item{\code{Age}}{a numeric vector}
    \item{\code{Sex}}{a factor with levels \code{0} \code{1}}
    \item{\code{ChestPainType}}{a factor with levels \code{1} \code{2} \code{3} \code{4}}
    \item{\code{RestBloodPressure}}{a numeric vector}
    \item{\code{SerumCholestoral}}{a numeric vector}
    \item{\code{FastingBloodSugar}}{a factor with levels \code{0} \code{1}}
    \item{\code{ResElectrocardiographic}}{a factor with levels \code{0} \code{1} \code{2}}
    \item{\code{MaxHeartRate}}{a numeric vector}
    \item{\code{ExerciseInduced}}{a factor with levels \code{0} \code{1}}
    \item{\code{Slope}}{a factor with levels \code{1} \code{2} \code{3}}
    \item{\code{MajorVessels}}{a factor with levels \code{0} \code{1} \code{2} \code{3}}
    \item{\code{Thal}}{a factor with levels \code{3} \code{6} \code{7}}
    \item{\code{Class}}{a factor with levels \code{1} \code{2}}
  }
}

\examples{
data(heart)
}
\keyword{datasets}
