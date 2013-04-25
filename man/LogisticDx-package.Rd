\name{LogisticDx-package}
\alias{LogisticDx-package}
\alias{LogisticDx}
\alias{logisticDx}
\docType{package}
\title{
  Diagnostic tests for logistic regression models
}
\description{
  Generate sample data frames for logistic regression. Diagnostic tests and plots for logistic regression models
}
\details{
  \tabular{ll}{
    Package: \tab LogisticDx\cr
    Type: \tab Package\cr
    Version: \tab 0.5\cr
    Date: \tab 2012-09-03\cr
    Date: \tab 2012-08-26 \cr
    License: \tab GPL (>= 2) \cr
    LazyLoad: \tab yes \cr
  }
  \item {\code{genLogiDf}}{generate data frame for logistic regression and model for the data frame}
  \item {\code{logiDx}}{return matrix of diagnostics for a logistic model}
  \item {\code{logiFactLev}}{return contingency tables for predictors in model which are factrs; also plots these}
  \item {\code{logiGOF}}{return goodness of fit restults for a logistic model}
  \item \code{logiProb} {return odds ratios and probabilities for all combinations of predictors in model}
  \item {\code{logiSS}}{return \dfn{n} the number of observations required to verify that the value of a predicted coefficient is that determined by the model, with given power and sensitivity}
  \item {\code{plotLogiDx}}{ diagnostic plots illustrating key diagnostics for a logistic model}
  \item {\code{stukel}}{Stukel's diagnostic test for the logistic fit}

  }
\author{
  Christopher Dardis <christopherdardis@gmail.com>
}
\references{
  Hosmer, D. and Lemeshow, S 2000 \emph{Applied Logistic Regression. Second edition}. John Wiley and Sons, Inc.
}
\keyword{package}
\examples{
model1 <- genLogiDf(n=20)$model
logiDx(model1)
}
