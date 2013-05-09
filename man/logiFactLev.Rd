\name{logiFactLev}
\alias{logiFactLev}
\title{Plots factors and gives tables of factors by level}
\usage{
  logiFactLev(df, lev = 4, alpha = 0.05, y = "y",
    tukey = FALSE, round1 = 2, width = 1366, height = 768)
}
\arguments{
  \item{df}{A dataframe, with predictor variables (binary,
  factors and continuous) and a binary outcome}

  \item{lev}{no. of \dfn{levels} which defines a factor.
  (I.e. column in data frame is a factor if it has 2<
  levels \eqn{\leq}{<=} \emph{lev})}

  \item{alpha}{Significance level \eqn{\alpha}{alpha}}

  \item{y}{Name of the outcome variable}

  \item{tukey}{If \code{TRUE}, will add letters below plots
  giving Tukeys all-pairwise comparison; identical letters
  indicate means are likely to be similar by this test (at
  p<0.05)}

  \item{round1}{No. digits to which to round (for display)}

  \item{width}{Width of screen(display device) in pixels}

  \item{height}{Width of screen(display device) in pixels}
}
\value{
  Returns plots as detailed above, one for each factor,
  each on a new screen. \cr \cr Also returns a list, with
  one level for each factor. Each item on the list is a
  matrix showing numbers of outcome \eqn{y=0} and \eqn{y=1}
  by level of factor as well as \% where \eqn{y=1}, with
  1-\code{alpha}\%CI given from Jeffreys prior.
}
\description{
  For a given data frame, will check which predictors are
  \dfn{factors} i.e. those with more with no. levels >2 and
  <\code{lev}. \cr \cr For each factor a plot is returned
  with \% where \emph{y}=1. A (1-\code{alpha})\% confidence
  interval is given using Jerrreys prior. \cr \cr A
  contingency table of outcome by factor level is plotted
  and returned. \cr \cr Tukeys test is generated from
  \code{multcomp::glht}.
}
\note{
  The letter display for Tukeys all-pairwise comparison
  relies on \code{multcomp::cld} which can be slow for
  factors with e.g. >20 levels. Plots work best for factors
  with small number of levels and can become
  \sQuote{overcrowded} with e.g. >20 levels.
}
\examples{
set.seed(1)
df1 <- genLogiDf(b=1,c=1,f=3,n=20)$df
logiFactLev(df1,tukey=TRUE)
}
\seealso{
  \code{\link{glht}} \code{\link{cld}}
}
\keyword{hplot}
\keyword{univar}

