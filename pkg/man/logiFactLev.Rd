\name{logiFactLev}
\alias{logiFactLev}
\title{
Plots factors and gives tables of factors by level
}
\description{
For a given data frame, will check which predictors are \dfn{factors}
i.e. those with more than two levels and the ratio of levels to number
of observations is < \emph{r}. (I.e. number of levels is < \emph{r *
  number observations}).
For each factor a plot is returned showing the level of the factor, with
\% where outcome \emph{y}=1. 95\% confidence interval given from
Jerrreys prior. A contingency table of outcome by factor level is
plotted and returned.
}
\usage{
logiFactLev(df, r = 0.4, y = "y",
tukey = FALSE, width = 1366, height = 768)
}
\arguments{
  \item{df}{
A dataframe, with predictor variables (binary, factors and continuous)
and a binary outcome.
}
  \item{r}{
The ratio of levels to the number of observations (number of rows in
data frame). (I.e. number of levels is < \emph{r * number
  observations}).
}
  \item{y}{
Name of the outcome variable.
}
  \item{tukey}{
If \code{TRUE}, will add letters below plots giving Tukeys all-pairwise
comparison; identical letters indicate means are likely to be similar by
this test (at p<0.05).
}
  \item{width}{
Width of display device/ screen (pixels).
}
  \item{height}{
Width of display device/ screen (pixels).
}
}
\details{
Tukeys test is generated from \pkg{multcomp::glht}.
}
\value{
  Returns plots as detailed above, one for each factor, each on a new screen.
  A list, with one level for each factor. Each item on the list is a
  matrix showing numbers of outcome \eqn{y=0} and \eqn{y=1} by level of factor as well
  as \% where \eqn{y=1}, with 95\%CI given from Jeffreys prior.
}
\note{
The letter display for Tukeys all-pairwise comparison relies on
\pkg{multcomp::cld} which can be slow for factors with e.g. >20 levels.
Plots work best for factors with small number of levels and can become
\sQuote{overcrowded} with e.g. >20 levels.
}
\seealso{
  \code{\link{glht}}
  \code{\link{cld}}
}
\examples{
\dontrun{
df1 <- genLogiDf(b=1,c=1,f=3,n=20)$df
logiFactLev(df1)
  }
}
\keyword{hplot}
\keyword{univar}
