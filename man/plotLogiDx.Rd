\name{plotLogiDx}
\alias{plotLogiDx}
\title{Diagnostic plots for a logistic regression}
\usage{
  plotLogiDx(model, noPerPage = 12, cexp = 2,
    cols = c("deepskyblue", "dodgerblue"),
    identify = FALSE, extras = FALSE, width = 1500,
    height = 800)
}
\arguments{
  \item{model}{A logistic regression model of class
  \bold{glm}}

  \item{noPerPage}{Number of plots per page (for initial
  plots). Will be used as \emph{guidance} and optimised for
  ease of display}

  \item{cols}{Colours. Used by \code{graphics::points}}

  \item{cexp}{Cex (Character EXpansion). Used by
  \code{graphics::points}}

  \item{identify}{If \code{TRUE} will give option to
  identify individual points on a number of the plots
  produced.}

  \item{extras}{If \code{TRUE} produces additional plots,
  detailed below}

  \item{width}{Width of screen(display device) in pixels}

  \item{height}{Height of screen(display device) in pixels}
}
\value{
  The following are plotted, for each covariate group:

  \item{probabiility_x_leverage}{Probability of \emph{y=1}
  for this group by leverage (diagonal of hat matrix, a
  measure of influence)}
  \item{probability_x_dXsq}{\dfn{dXsq} change in Pearson
  chi-square statistic with deletion of this group}
  \item{probability_x_dBhat}{\dfn{dBhat} change in Bhat;
  the difference in the maximum likelihood estimators
  \bold{Beta} for model coefficients with all subjects
  included vs those with this group, standardized by the
  estimated covariance matrix of \bold{Beta}}
  \item{probability_x_dDev}{\dfn{dDev} change in Deviance
  when this group excluded} \item{bubbleplot}{Probability
  by dXsq, with Area proportional to dBhat}
  \item{leverage_x_dXsq}{\dfn{dXsq} change in Pearson
  chi-square statistic with deletion of this group}
  \item{leverage_x_dBhat}{\dfn{dBhat} change in Bhat; the
  difference in the maximum likelihood estimators
  \bold{Beta} for model coefficients with all subjects
  included vs those with this group, standardized by the
  estimated covariance matrix of \bold{Beta}}
  \item{leverage_x_dDev}{\dfn{dDev} change in Deviance when
  this group excluded} \item{ROC}{Receiver Operator Curve}

  Additional plots given by \code{extras=TRUE}:

  \item{influenceplot}{from \code{\link{influencePlot}}}
  \item{studentizedResiduals_x_hatvalues}{\dfn{Studentized
  residual} = residual / estimate of its standard
  deviation} \item{spreadlevelplot}{from
  \code{\link{spreadLevelPlot}}}
  \item{qqPlot}{quantile-quantile plot vs Normal for
  residuals} \item{influenceIndexPlot}{Cooks distance,
  studentized residual and hat values for each observation}
  \item{pairsplot}{for measure of influence dBhat, dXsq,
  dDev} \item{component+residualplots}{from
  \code{\link{cr.plot}}} \item{added-variableplots}{from
  \code{\link{av.plot}}} \item{marginalmodelplots}{from
  \code{\link{marginalModelPlot}}}
}
\description{
  Common diagnostic plots for a logistic regression model
}
\note{
  Different colors can be found with e.g.
  \code{grDevices::colours()[grep("blue",grDevices::colours())]}
}
\examples{
set.seed(1)
### generate 8x covariate patterns

mod1 <- genLogiDf(b=3,f=0,c=0,n=50)$model
plotLogiDx(mod1, cexp=8, noPerPage=1)
plotLogiDx(mod1, cexp=3, noPerPage=6, extras=TRUE)
}
\seealso{
  \code{\link{car}}
}
\keyword{hplot}

