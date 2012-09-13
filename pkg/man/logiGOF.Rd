\name{logiGOF}
\alias{logiGOF}
\title{
Goodness of fit tests for a logistic regression model
}
\description{
Gives 15 commonly employed measures of goodness of fit for a logistic
regression model
}
\usage{
logiGOF(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{A model of class \bold{glm}}
}

\value{
  \item{ChiPearCov}{Pearsons Chi-square, calculated by covariate group,
    with p value and interpretation}
  \item{ChiPearIndiv}{Pearsons Chi-square, calculated by individual observation,
    with p value and interpretation}
  \item{ChiPearTab}{Pearsons Chi-square, calculated by table of
    covariate patterns by outcome, with p value and interpretation}
  \item{OsRo}{Osius & Rojek test of the logistic link, with p value and interpretation}
  \item{ChiDevCov}{Deviance Chi-square, calculated by covariate group,
    with p value and interpretation}
  \item{ChiDevIndiv}{Deviance Chi-square, calculated by individual observation,
    with p value and interpretation}
  \item{ChiDevTab}{ Deviance Chi-square, calculated by table of
    covariate patterns by outcome, with p value and interpretation}
  \item{CovPatTab}{ Matrix of covariance patterns, used to
      calculate above Chi-square tests of Pearson residuals and Deviance}
  \item{HosLem}{Hosmer & Lemeshow goodness of fit test, with \emph{g=10} quantile groups,with p value and interpretation}
  \item{modHosLem}{modified Hosmer & Lemeshow goodness of fit test, with \emph{g=10} quantile groups, with p value and interpretation}
  \item{CesHou}{le Cessie, van Houwelingen, Copas & Hosmer unweighted
	sum of squares test for global goodness of fit, with p value and
	interpretation}
  \item{Stuk}{Stukels test of the appropriateness of the logistic
	link, with p value and interpretation}
      \item{PR2}{Pearsons R^2, correlation of observed outcome with
	predicted}
      \item{ssR2}{linear regression-like sum of squares R^2, using
	covariate patterns}
      \item{llR2}{log-likelohood based R^2, calculated by covariate group}
      \item{ROC}{Area under the Receiver Operating Curve, with 95\% CI by
	method of DeLong}
    }
\author{
Yongmei Ni: modified Hosmer & Lemeshow goodness of fit test.
}
\section{Warning}{
  Will fail if cannot generate a hat matrix for the model using \code{logiDx.R}
}
\seealso{
\code{\link{logiDx}}
}
\examples{
\dontrun{
f1 <- genLogiDf(n=100)$model
logiGOF(f1)
  }
}
\keyword{htest}