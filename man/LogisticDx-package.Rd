\docType{package}
\name{LogisticDx-package}
\alias{logisticDx}
\alias{LogisticDx}
\alias{LogisticDx-package}
\title{Diagnostic tests for logistic regression models}
\description{
  Generate sample data frames for logistic regression.
  Diagnostic tests and plots for logistic regression models
}
\details{
  \tabular{ll}{ Package: \tab LogisticDx\cr Type: \tab
  Package\cr Version: \tab 0.5\cr Date: \tab 2012-09-03\cr
  Date: \tab 2012-08-26 \cr License: \tab GPL (>= 2) \cr
  LazyLoad: \tab yes \cr }
}
\examples{
model1 <- genLogiDf(n=20)$model
logiDx(model1)
}
\author{
  Christopher Dardis <christopherdardis@gmail.com>
}
\references{
  Hosmer, D. and Lemeshow, S 2000 \emph{Applied Logistic
  Regression. Second edition}. John Wiley and Sons, Inc.
}
\keyword{package}

