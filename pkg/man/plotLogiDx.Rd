\name{plotLogiDx}
\alias{plotLogiDx}
\title{
Diagnostic plots for a logistic regression
}
\description{
Common diagnostic models for a logistic regression model.
}
\usage{
plotLogiDx( model, noPerPage = 12, identify = FALSE,
extras = FALSE, width = 1500, height = 800)
}
\arguments{
  \item{model}{
    A logistic regression model of class \bold{glm}
}
  \item{noPerPage}{
    Number of plots per page (for initial plots)
}
  \item{identify}{
    If \code{identify=TRUE} will give option to identify individual
    points on a number of the plots produced.
}
  \item{extras}{
    If \code{extras=TRUE} produces additional plots
    }
  \item{width}{
    Width of display device/ screen in pixels
}
  \item{height}{
    Height of display device/ screen in pixels
}
}
\value{
  The following are plotted, for each covariate group:
  \item{probabiility by leverage}{Probability of \emph{y=1} for this group by leverage (diagonal of hat matrix, a measure of influence)}
  \item{probability by dXsq}{\dfn{dXsq} change in Pearson chi-square statistic with deletion of this group}
  \item{probability by dBhat}{\dfn{dBhat} change in Bhat; the difference in the maximum likelihood estimators \bold{Beta} for model coefficients with all subjects included vs those with this group, standardized by the estimated covariance matrix of \bold{Beta}}
  \item{probability by dDev}{\dfn{dDev} change in Deviance when this group excluded}
  \item{bubble plot}{Probability by dXsq, with Area proportional to dBhat}
  \item{leverage by dXsq}{\dfn{dXsq} change in Pearson chi-square statistic with deletion of this group}
  \item{leverage by dBhat}{\dfn{dBhat} change in Bhat; the difference in
    the maximum likelihood estimators \bold{Beta} for model coefficients with all subjects included vs those with this
    group, standardized by the estimated covariance matrix of
    \bold{Beta} }
  \item{leverage by dDev}{\dfn{dDev} change in Deviance when this group excluded}
  \item{ROC}{Receiver Operator Curve}
  \item{influence plot}{from \code{\link{influencePlot}}}

  Additional plots given by \code{extras=TRUE}:

  \item{studentized residuals by hat values}{Studentized residual = residual / estimate of its standard deviation}
  \item{spread level plot}{from \code{\link{spreadLevelPlot}}}
  \item{qqPlot}{quantile-quantile plot vs Normal for residuals}
  \item{influenceIndexPlot}{Cooks distance, studentized residual and hat values for each observation}
  \item{pairs plot}{for measure of influence dBhat, dXsq, dDev}
  \item{component + residual plots}{from \code{\link{cr.plot}}}
  \item{added-variable plots}{from \code{\link{av.plot}}}
  \item{marginal model plots}{from \code{\link{marginalModelPlot}}}
  }
\seealso{
\code{\link{car}}
}

\examples{
\dontrun{
f1 <- genLogiDf(b=3,f=0,c=0,n=50)$model
plotLogiDx(f1) # 8 covariate patterns
   }
}
\keyword{hplot}

