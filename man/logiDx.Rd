\name{logiDx}
\alias{logiDx}
\title{Diagnostics for logistic regression}
\usage{
  logiDx(model, width = TRUE, round1 = 3)
}
\arguments{
  \item{model}{A model of class \bold{glm}}

  \item{width}{if \code{width = TRUE} changes
  \code{options(width)} to allow all columns in result to
  be displayed side-by-side}

  \item{round1}{No. digits to which to round (for display)}
}
\value{
  A list containing one matrix, \dQuote{dxMatrix}. \cr \cr
  The initial columns show all unique combinations of the
  predictor variables (i.e. one row for each
  \emph{covariate pattern}). \cr \cr Subsequent columns are
  labelled as follows:

  \item{obs}{Number of observations with this covariate
  pattern} \item{prob}{Probability of this covariate
  pattern} \item{yhat}{Number of observations of \eqn{y=1},
  predicted by the model} \item{y}{\emph{Actual} number of
  observations of \eqn{y=1} from the data}
  \item{lev}{\dfn{leverage}, the diagonal of the hat matrix
  used to generate the model; a measure of influence of
  this covariate pattern} \item{devR}{\dfn{Deviance
  residual}, calculated by covariate pattern; a measure of
  influence of this covariate pattern}
  \item{PeR}{\dfn{Pearson residual}, calculated by
  covariate pattern; a measure of influence of this
  covariate pattern. Given by: \deqn{
  \sqrt{obs}\sqrt{\frac{prob}{(1-prob)}}}{ obs^0.5
  (prob/1-prob)^0.5}} \item{sPeR}{\dfn{standardized Pearson
  residual} calculated by covariate pattern; a measure of
  influence of this covariate pattern. Given by: \deqn{
  \frac{PeR}{\sqrt{(1-lev)}}}{ PeR.(1-lev)^0.5}}
  \item{dBhat}{\dfn{change in Bhat}, the standardized
  difference between the original maximum likelihood
  estimates \bold{B} and that the estimates with this
  covariate pattern excluded} \item{dXsq}{\dfn{change in
  Chi-square}, decrease in the value of Pearson chi-square
  statistic with this covariate pattern excluded. Given by:
  \deqn{sPeR^2}} \item{dDev}{\dfn{change in deviance}
  \bold{D} with this covariate pattern excluded. Given by:
  \deqn{ \frac{dev^2}{(1-lev)}}{ d^2/(1-lev)}}
}
\description{
  Returns standard diagnostic measures for a logistic
  regression model
}
\note{
  Values for the statistics are calculated by
  \emph{covariate pattern}. Different values may be
  obtained if calculated for each individual obervation
  (i.e. row in data frame). \cr \cr Generally, the values
  calculated by covariate pattern are preferred,
  particularly where \eqn{obs >5}.
}
\examples{
m1 <- genLogiDf()$model
logiDx(m1)
}
\seealso{
  \code{\link{plotLogiDx}}
}
\keyword{array}

