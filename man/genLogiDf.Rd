\name{genLogiDf}
\alias{genLogiDf}
\title{Generate logistic data frame and model}
\usage{
  genLogiDf(f = 0, b = 2, c = 1, n = 20, rb = 0.5, nf = 3,
    rc = 0.8, ry = 0.5, asFactor = FALSE, timelim = 5,
    speedglm = FALSE)
}
\arguments{
  \item{b}{\dfn{binomial predictors}, the number of
  predictors which are binary, i.e. limited to 0 or 1}

  \item{f}{\dfn{factors}, the number of predictors which
  are factors}

  \item{c}{\dfn{continuous predictors}, the number of
  predictors which are continuous}

  \item{n}{number of observations in the data frame}

  \item{nf}{the no. of levels in a factor}

  \item{rb}{\dfn{ratio for binomnial predictors} the ratio
  of 1s to total observations for the binomial predictors
  e.g. if \code{rb=0.3}, 30\% will be 1s, 70\% will be 0s}

  \item{rc}{\dfn{ratio for continuous variables} the ratio
  of levels of continuous variables to the total number of
  observations \dfn{n} e.g. if \code{rc=0.8} and
  \code{n=100}, it will be in the range 1-80}

  \item{ry}{\dfn{ratio for y} the ratio of 1s to total
  observations for the binomial predictors e.g. if
  \code{ry=0.5}, 50\% will be 1s, 50\% will be 0s}

  \item{asFactor}{if \code{TRUE}, predictors given as
  factors will be converted to factors in data frame before
  model is fit}

  \item{timelim}{function will timeout after \code{timelim}
  secs}

  \item{speedglm}{return fitted model with \code{speedglm}
  instead of \code{glm}}
}
\value{
  A list with the following values: \item{df}{data frame
  with predictors (labelled \eqn{x1,x2, ..., xn}) and
  outcome (y), with \emph{n} rows (observations)}
  \item{model}{model fit with \code{
  glm(family=binomial())} or
  \code{speedglm(family=binomial())} }
}
\description{
  Generates a data frame with a binary outcome, and a
  logistic model to describe it. Model is fitted with
  \code{glm()}.
}
\note{
  Using \code{asFactor=TRUE} with factors which have a
  large number of levels (e.g. \code{nf >30}) on large
  datasets (e.g. \eqn{n >1000}) can cause fitting to be
  excessively slow.
}
\examples{
set.seed(1)
genLogiDf()
genLogiDf(b=0,c=2,n=100,rc=0.7)
genLogiDf(b=1,c=0,n=1000)
}
\keyword{datagen}

