\name{logiGOF}
\alias{logiGOF}
\title{Goodness of fit tests for a logistic regression model}
\usage{
  logiGOF(model, g = 10)
}
\arguments{
  \item{model}{A model of class \bold{glm}}

  \item{g}{No. groups (quantiles) into which to split
  observatiosn for Hosmer-Lemeshow and modified
  Hosmer-Lemeshow tests.}
}
\value{
  A list with the following items:
  \item{ChiPearCov}{Pearsons Chi-square, calculated by
  covariate group, with p value and interpretation}
  \item{ChiPearIndiv}{Pearsons Chi-square, calculated by
  individual observation, with p value and interpretation}
  \item{ChiPearTab}{Pearsons Chi-square, calculated by
  table of covariate patterns by outcome, with p value and
  interpretation} \item{OsRo}{Osius & Rojek test of the
  logistic link, with p value and interpretation}
  \item{ChiDevCov}{Deviance Chi-square, calculated by
  covariate group, with p value and interpretation}
  \item{ChiDevIndiv}{Deviance Chi-square, calculated by
  individual observation, with p value and interpretation}
  \item{ChiDevTab}{Deviance Chi-square, calculated by table
  of covariate patterns by outcome, with p value and
  interpretation} \item{CovPatTab}{Matrix of covariance
  patterns, used to calculate above Chi-square tests of
  Pearson residuals and Deviance} \item{HosLem}{Hosmer &
  Lemeshow goodness of fit test, with \emph{g=10} quantile
  groups,with p value and interpretation}
  \item{modHosLem}{modified Hosmer & Lemeshow goodness of
  fit test, with \emph{g=10} quantile groups, with p value
  and interpretation} \item{CesHou}{le Cessie, van
  Houwelingen, Copas & Hosmer unweighted sum of squares
  test for global goodness of fit, with p value and
  interpretation} \item{Stuk}{Stukels test of the
  appropriateness of the logistic link, with p value and
  interpretation} \item{PR2}{Pearsons R^2, correlation of
  observed outcome with predicted} \item{ssR2}{linear
  regression-like sum of squares R^2, using covariate
  patterns} \item{llR2}{log-likelohood based R^2,
  calculated by covariate group} \item{ROC}{Area under the
  Receiver Operating Curve, with 95\% CI by method of
  DeLong}
}
\description{
  Gives 15 commonly employed measures of goodness of fit
  for a logistic regression model
}
\note{
  Warning: Will fail if cannot generate a hat matrix for
  the model using \code{logiDx}
}
\examples{
set.seed(1)
m1 <- genLogiDf(n=100)$model
}
\author{
  Yongmei Ni: modified Hosmer & Lemeshow goodness of fit
  test (adapted)
}
\seealso{
  \code{\link{logiDx}}
}
\keyword{htest}

