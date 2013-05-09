\name{logiProb}
\alias{logiProb}
\title{Logits, odds ratios and probabilities for all
combinations of coefficients}
\usage{
  logiProb(model)
}
\arguments{
  \item{model}{A logistic regression model of class
  \bold{glm}}
}
\value{
  A matrix giving, for each combination of predictors
  (coefficients): \item{logit}{The logit for a given
  combination of coefficients} \item{ln(OR)}{Natural log of
  odds ratio} \item{OR}{Odds ratio} \item{p}{probability}
}
\description{
  For all combinations of predictors in model (including
  intercept only), generate logit, log odds ratio (OR), OR
  and probability.

  Values are calculated for a change in the value of the
  coeffient for the predictor from 0 to 1. (For continuous
  predictors changes of more than one unit may have more
  practical significance).
}
\examples{
set.seed(1)
f1 <- genLogiDf(n=50)$model
logiProb(f1)
}
\keyword{array}
\keyword{math}

