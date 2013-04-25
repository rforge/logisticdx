\name{stukel}
\alias{stukel}
\title{
Stukels test of the logistic link
}
\description{
Calculates Stukels test for a logistic regression model. This determines
the appropriateness of the logistic link.
Two new covariates, z1 and z2 are generates such that
z1=0.5*logit^2*I(pi>=0.5), z2=-0.5*logit^2*I(pi<=0.5) where I(arg) is 1
where arg is true and 0 where false. Tests if significant change in the
model with the addition of these coefficients.
}
\usage{
stukel(object, alternative = c("both", "alpha1", "alpha2"))
}
\arguments{
  \item{object}{
    An object of class \bold{glm}
  }
  \item{alternative}{add both z1 and z2 to model or just one of the above}
   }
\value{
  \item{statistic}{value of statistic}
  \item{parameter}{degrees of freedom}
  \item{p.value}{if <0.05 suggests should accept null hypothesis that logistic model is incorrect for the data}
  \item{alternative}{alternative}
  \item{method}{method}
  \item{data.name}{name of object}
  \item{allstat}{statistics for all tests}
  \item{allpar}{degrees of freedom}
  \item{allpval}{all p values}
 }
\references{
\href{http://www.stat.ufl.edu/~presnell/Courses/sta7249-2008sp/R/Src/}{source code}
}
\author{
Brett Presnell
}
\examples{
\dontrun{
f1 <- genLogiDf(b=3,f=0,c=0,n=50)$model
stukel(f1)
  }
}
\keyword{htest}