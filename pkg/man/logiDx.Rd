\name{logiDx}
\alias{logiDx}
\title{
  Diagnostics for logistic regression
}
\description{
  Returns standard diagnostic measures for a logistic regression model
}
\usage{
logiDx(model, width=TRUE)
}
\arguments{
  \item{model}{A model of class \bold{glm}}
  \item{width}{if \code{width = TRUE} changes \code{options(width)} to allow all columns in result to be displayed side-by-side}
}
\value{
  A list containing one matrix, \dQuote{dxMatrix}.
  The initial columns show all combinations of the predictor variables
  (i.e. one for each \emph{covariate pattern}). Subsequent columns are
  labelled as follows:
  \item{obs}{Number of observations with this covariate pattern}
  \item{pi}{Probability \eqn{\pi} of this covariate pattern}
  \item{yhat}{Number of observations of \eqn{y=1}, predicted by the model}
  \item{y}{\emph{Actual} number of observations of \eqn{y=1} from the data}
  \item{lev}{\dfn{leverage}, the diagonal of the hat matrix used to generate the model; a measure of influence of this covariate pattern}
  \item{devR}{\dfn{Deviance residual}, calculated by covariate pattern; a measure of influence of this covariate pattern}
  \item{PeR}{\dfn{Pearson residual}, calculated by covariate pattern; a measure of influence of this covariate pattern. Given by: \deqn{ \sqrt{obs}\sqrt{\frac{\pi}{(1-\pi)}}} }
  \item{sPeR}{\dfn{standardized Pearson residual} calculated by covariate pattern; a measure of influence of this covariate pattern. Given by: \deqn{ \frac{PeR}{\sqrt{(1-lev)}}} }
  \item{dBhat}{\dfn{change in Bhat}, the standardized difference between the original maximum likelihood estimates \bold{B} and that the estimates with this covariate pattern excluded}
  \item{dXsq}{\dfn{change in Chi-square}, decrease in the value of Pearson chi-square statistic with this covariate pattern excluded. Given by: \deqn{sPeR^2} }
  \item{dDev}{\dfn{change in deviance} \bold{D} with this covariate pattern excluded. Given by: \deqn{ \frac{d^2}{(1-lev)}}}
  }
\note{
Values for the statistics are calculated by \emph{covariate
  pattern}. Different values may be obtained if calculated for each
individual obervation (i.e. row in data frame). Generally, the values calculate by covariate pattern are preferred, particularly where
\dfn{obs} is >5.
}
\seealso{
\code{\link{plotLogiDx}}
}
\examples{
m1 <- genLogiDf()$model
logiDx(m1)
}
\keyword{array}
