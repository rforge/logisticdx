\name{genLogiDf}
\alias{genLogiDf}
\title{
Generate logistic data frame and model
}
\description{
Generates a data frame with a binary outcome, and a logistic model to
describe it. Model if fitted with \code{glm}
}
\usage{
genLogiDf (b = 2, f = 0, c = 1, n = 10, rb = 0.5,
rf = 0.3, rc = 0.8, ry = 0.5,
asFactor = FALSE, timeout = 5)
}
\arguments{
  \item{b}{\dfn{binomial predictors}, the number of predictors which are binary, i.e. limited to 0,1 in range}
  \item{f}{\dfn{factors}, the number of predictors which are factors}
  \item{c}{\dfn{continuous predictors}, the number of predictors which are continuous}
  \item{n}{number of observations in the data frame}
  \item{rb}{\dfn{ratio for binomnial predictors} the ratio of 1s to total observations for the binomial predictors e.g. if = 0.3, 30\% will be 1s, 70\% will be 0s}
  \item{rf}{\dfn{ratio for factors} the ratio of levels of factor to the total number of observations \dfn{n} e.g. if = \code{rf = 0.3} and \code{n=100}, factor will have potentially 30 levels, i.e. 1-30}
  \item{rc}{\dfn{ratio for continuous variables} the ratio of levels of
    continuous variables  to the total number of observations \dfn{n}
    e.g. if \code{rc = 0.8} and \code{n}=100, factor will have potentially 80 levels, i.e. 1-80}
  \item{ry}{\dfn{ratio for y}  the ratio of 1s to total observations for the binomial predictors e.g. if \code{ry = 0.5}, 50\% will be 1s, 50\% will be 0s}
  \item{asFactor}{if \code{TRUE}, predictors given as factors will be converted to factors in data frame before model is fit}
  \item{timeout}{function will timeout after this variable elapses}
  }

\value{
  A list with the following values:
  \item{df}{ data frame with predictors (labelled x1,x2...xn) and outcome
    (y), with \emph{n} observations}
  \item{model}{ model fit with with \code{ glm(family = binomial())}
    }
}
\note{
Using \code{asFactor=TRUE} with factors which have a large number of
levels (e.g. \code{rf>0.5}) on large datasets (e.g. \code{n >1000}) can cause
fitting to be excessively slow.
}
\examples{
genLogiDf()
genLogiDf(b=0,c=2,n=100,rc=0.7)
}
\keyword{datagen}
