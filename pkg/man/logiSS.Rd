\name{logiSS}
\alias{logiSS}
\title{
Sample size for given coefficient and events per covariate for model
}
\description{
Gives sample size necessary to demonstrate that coefficient in model for given
predictor is equal to its given value (rather than equal to zero) for a given level
of power and significance. Also number of events (smaller of outcome \eqn{y=0}
and outcome \eqn{y=1}) per predictor.
}
\usage{
logiSS( model, alpha = 0.05, beta = 0.8, coeff = "x1")
}
\arguments{
  \item{model}{
   A logistic regression model of class \bold{glm}
    }
  \item{alpha}{
    significance level \eqn{\alpha} for null-hypothesis significance test
  }
  \item{beta}{
     power \eqn{\beta} for null-hypothesis significance test
}
  \item{coeff}{
    Name of predictor (coefficient) in model to be tested
  }
}
\details{
Uses different methods depending on whether model has one binomial, one
continuous or multiple predictors.
}
\value{
  A list of:
    \item{res}{Result: Sample size required to show coefficient for
      predictor is as given in the model rather than 0}
    \item{epc}{Events per covariate; should be >10 to make meaningful
      statements about coefficients obtained}

}
\examples{
f1 <- genLogiDf(n=50)$model
logiSS(f1)
}
\keyword{htest}
