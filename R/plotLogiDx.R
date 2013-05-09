##' @name plotLogiDx
##' @export
##' @include logiDx.R
##' @title Diagnostic plots for a logistic regression
##'
##' @description
##' Common diagnostic plots for a logistic regression model
##'
##' @param model A logistic regression model of class \bold{glm}
##' @param noPerPage Number of plots per page (for initial plots).
##' Will be used as \emph{guidance} and optimised for ease of display
##' @param cols Colours. Used by \code{graphics::points}
##' @param cexp Cex (Character EXpansion). Used by \code{graphics::points}
##' @param identify If \code{TRUE} will give option to identify
##' individual points on a number of the plots produced.
##' @param extras If \code{TRUE} produces additional plots, detailed below
##' @param width Width of screen(display device) in pixels
##' @param height Height of screen(display device) in pixels
##' @return The following are plotted, for each covariate group:
##'
##' \item{probabiility_x_leverage}{Probability of \emph{y=1} for this group
##' by leverage (diagonal of hat matrix, a measure of influence)}
##' \item{probability_x_dXsq}{\dfn{dXsq} change in Pearson chi-square
##' statistic with deletion of this group}
##' \item{probability_x_dBhat}{\dfn{dBhat} change in Bhat; the difference in
##' the maximum likelihood estimators \bold{Beta} for model coefficients with
##' all subjects included vs those with this group, standardized by the
##' estimated covariance matrix of \bold{Beta}}
##' \item{probability_x_dDev}{\dfn{dDev} change in Deviance when this group
##' excluded}
##' \item{bubbleplot}{Probability by dXsq, with Area proportional to dBhat}
##' \item{leverage_x_dXsq}{\dfn{dXsq} change in Pearson chi-square statistic
##' with deletion of this group}
##' \item{leverage_x_dBhat}{\dfn{dBhat} change in Bhat; the difference in
##' the maximum likelihood estimators \bold{Beta} for model coefficients with
##' all subjects included vs those with this group, standardized by the
##' estimated covariance matrix of \bold{Beta}}
##' \item{leverage_x_dDev}{\dfn{dDev} change in Deviance when this group
##' excluded}
##' \item{ROC}{Receiver Operator Curve}
##'
##' Additional plots given by \code{extras=TRUE}:
##'
##' \item{influenceplot}{from \code{\link{influencePlot}}}
##' \item{studentizedResiduals_x_hatvalues}{\dfn{Studentized residual}
##' = residual / estimate of its standard deviation}
##' \item{spreadlevelplot}{from \code{\link{spreadLevelPlot}}}
##' \item{qqPlot}{quantile-quantile plot vs Normal for residuals}
##' \item{influenceIndexPlot}{Cooks distance, studentized residual and hat
##' values for each observation}
##' \item{pairsplot}{for measure of influence dBhat, dXsq, dDev}
##' \item{component+residualplots}{from \code{\link{cr.plot}}}
##' \item{added-variableplots}{from \code{\link{av.plot}}}
##' \item{marginalmodelplots}{from \code{\link{marginalModelPlot}}}
##'
##' @note Different colors can be found with e.g.
##' \code{grDevices::colours()[grep("blue",grDevices::colours())]}
##'
##' @seealso \code{\link{car}}
##' @keywords hplot
##' @examples
##'
##' set.seed(1)
##' ### generate 8x covariate patterns
##'
##' mod1 <- genLogiDf(b=3,f=0,c=0,n=50)$model
##' plotLogiDx(mod1, cexp=8, noPerPage=1)
##' plotLogiDx(mod1, cexp=3, noPerPage=6, extras=TRUE)
##'
plotLogiDx <- function(model, noPerPage=12,
                       cexp=2,
                       cols=c("deepskyblue","dodgerblue"),
                       identify=FALSE,
                       extras=FALSE, width=1500, height=800) {
    if (class(model)[1] !="glm")
        stop("plotLogiDx only applies to objects of class glm")
### get diagnostics for model
    m1 <- logiDx(model)$dxMatrix
### find best balance for noPerPage
    if (noPerPage==1) {
        nrow1 <- ncol1 <- 1
        } else {
            balance <- function(x) (x^2 + noPerPage )/ x
            nrow1 <- round(stats::optimize(balance,
                                           interval=seq(1:noPerPage))$minimum,
                           0)
            ncol1 <- round( noPerPage/nrow1, 0)
        }
### open plot window
    windows(record=TRUE, width=width, height=height)
    p <- par
### oma=outer margins, mar=margins, bottom,left,top,right
    par( mfrow=c(nrow1,ncol1), oma=c(0,0,4,0), mar=c(4,6,3,0.5) )
###
    maintext <- function(){
        tex1 <- paste0("Diagnostic plots for logistic regression \n ",
                       deparse(model$formula))
        graphics::mtext(tex1, line = 0.3, outer = TRUE)
        }
###-------------------------------------------
### probabiility by leverage
    graphics::plot(m1[ ,"prob"], m1[ ,"lev"],
                   xlab="Probability for this covariate pattern",
                   main="Probability by leverage",
                   ylab="Leverage (hat matrix diagonal)")
### to make horizontal on y-axis:
### graphics::mtext("Leverage \n (hat \n matrix \n diagonal)", side=2, line=3, las=1)
    graphics::points( m1[ ,"prob"], m1[ ,"lev"], pch=21, cex=cexp,
           col = cols, bg = cols)
    if (isTRUE(identify)) graphics::identify(m1[ ,"prob"], m1[ ,"lev"])
    maintext()
###-------------------------------------------
### probability by ...
### dXsq
    graphics::plot(m1[ ,"prob"], m1[ ,"dXsq"],
                   xlab="Probability for this covariate pattern",
                   ylab="dXsq = decrease in Pearson Chi-sq \n without this pattern",
                   main="dXsq by probability")
    graphics::points(m1[ ,"prob"], m1[ ,"dXsq"], pch=21, cex=cexp,
           col = cols, bg=cols)
    if (isTRUE(identify)) graphics::identify(m1[ ,"prob"], m1[ ,"dXsq"])
    maintext()
### dBhat
    graphics::plot(m1[ ,"prob"], m1[ ,"dBhat"],
                   xlab="Probability for this covariate pattern",
                   ylab="dBhat = decrease in Bhat \n without this pattern",
                   main="dBhat by probability")
    graphics::points(m1[ ,"prob"], m1[ ,"dBhat"], pch=21, cex=cexp,
           col = cols, bg=cols)
    if (isTRUE(identify)) graphics::identify(m1[ ,"prob"], m1[ ,"dBhat"])
    maintext()
### dDev
    graphics::plot(m1[ ,"prob"], m1[ ,"dDev"],
                   xlab="Probability for this covariate pattern",
                   ylab="dDev = decrease in Deviance \n without this pattern",
                   main="dDev by probability")
    graphics::points(m1[ ,"prob"], m1[ ,"dDev"], pch=21, cex=cexp,
           col = cols, bg=cols)
    if (isTRUE(identify)) graphics::identify(m1[ ,"prob"], m1[ ,"dDev"])
    maintext()
###-------------------------------------------
### bubble plot - prob by dXsq, with area = dBhat
    radius <- sqrt(m1[ ,"dBhat"]/ m1[ ,"prob"])
    graphics::symbols(m1[ ,"prob"], m1[ ,"dXsq"],
                      circles=radius, inches=0.35,
                      fg="white", bg=cols,
                      xlab="Probability for this covariate pattern",
                      ylab="dXsq = decrease in Pearson Chi-sq \n without this pattern",
                      main = "Area proportional to dBhat \n Decrease in Bhat without this pattern")
    if (isTRUE(identify)) graphics::identify(m1[ ,"prob"], m1[ ,"dXsq"])
    maintext()

###---------------------------
### leverage by...
### dXsq
    graphics::plot(m1[ ,"lev"], m1[ ,"dXsq"],
                   xlab="Leverage for this covariate pattern",
                   ylab="dXsq = decrease in Pearson Chi-sq \n without this pattern",
                   main="dXsq by leverage")
    graphics::points(m1[ ,"lev"], m1[ ,"dXsq"], pch=21, cex=cexp,
           col = cols, bg=cols)
    if (isTRUE(identify)) graphics::identify(m1[ ,"lev"], m1[ ,"dXsq"])
    maintext()
### dBhat
    graphics::plot(m1[ ,"lev"], m1[ ,"dBhat"],
                   xlab="Leverage for this covariate pattern",
                   ylab="dBhat = decrease in Bhat \n without this pattern",
                   main="dBhat by leverage")
    graphics::points(m1[ ,"lev"], m1[ ,"dBhat"], pch=21, cex=cexp,
           col=cols, bg=cols)
    if (isTRUE(identify)) graphics::identify(m1[ ,"lev"], m1[ ,"dBhat"])
    maintext()
### dDev
    graphics::plot(m1[ ,"lev"], m1[ ,"dDev"],
                   xlab="Leverage for this covariate pattern",
                   ylab="dDev = decrease in Deviance \n without this pattern",
                   main="dDev by leverage")
    graphics::points(m1[ ,"lev"], m1[ ,"dDev"], pch=21, cex=cexp,
           col = cols, bg=cols)
    if (isTRUE(identify)) graphics::identify(m1[ ,"lev"], m1[ ,"dDev"])
    maintext()
###
### ROC curve
### get name of y/ outcome variable from formula
    r1 <- pROC::roc(response=model$data[[ncol(model$data)]],
                    predictor=model$fitted,
                     ci=TRUE, percent=TRUE)
### 0.5 = chance, aim >0.7
    pROC::plot.roc(r1, print.auc=TRUE, grid=TRUE,
                   print.auc.cex=0.8, main="ROC curve")
    maintext()
### influence plot
    if (isTRUE(identify)){
        car::influencePlot(model,
                           id.n=1.5,
                           id.method="identify",id.col="blue", scale=15,
                           xlab="Hat values (vertical lines at x2, x3 average hat value)",
                           main="Area proportional to Cooks distance",
                           sub="Cooks = change in coefficients if this point dropped")
    } else {
	car::influencePlot(model, id.n=1.5, id.method="noteworthy",
                           id.col="blue", scale=15,
                           xlab="Hat values (vertical lines at x2, x3 average hat value)",
                           main="Area proportional to Cooks distance",
                           sub="Cooks = change in coefficients if this point dropped")
    }
    maintext()
###
###---------------------------------------------------------
### optional plots
###
    if (isTRUE(extras)) {
### studentized residuals vs hat values
        graphics::plot(rstudent(model) ~ hatvalues(model),
                       xlab="Hat values",
                       ylab="Studentized residuals (residual / estimate of its S.D.)",
                       main="Studentized residuals by hat values")
        graphics::points(hatvalues(model), rstudent(model), pch=21, cex=cexp,
                         col = cols, bg=cols)
        if (isTRUE(identify)) {
            graphics::identify(hatvalues(model), rstudent(model), cex=1.5)}
        maintext()
### spread level plot
        car::spreadLevelPlot(model, layout=NA,
                             main="Spread-Level plot \n Spread should be even \n if variance of error is constant")
        maintext()
### qqnorm for residuals
        q1 <- stats::qqnorm(statmod::qresid(model),col="darkblue", cex=1.5,
                            main="Normal Q-Q plot for residuals from formula \n Check if residuals normally distributed",
                            ylab="Quantiles from formula")
        stats::qqline((statmod::qresid(model)),col="darkblue",lwd=1.5)
        if (isTRUE(identify)) graphics::identify(q1, cex=1.5)
        maintext()
### Cooks, stud~ res~ and hat values by index
        if (isTRUE(identify)){
            car::influenceIndexPlot(model, col="dodgerblue",
                                    id.method = "identify",
                                    vars=c("Cook", "Studentized", "hat"),
                                    main="Cooks distance, studentized residuals \n and hat values by Index (of observations)")
        } else {
            car::influenceIndexPlot(model, col="dodgerblue",
                                    id.method = "y",
                                    vars=c("Cook", "Studentized", "hat"),
                                    main="Cooks distance, studentized residuals \n and hat values by Index (of observations)")
        }
### pairs plot for measures of influence
        graphics::pairs(m1[ ,c("dBhat", "dXsq", "dDev")],
                        col="darkblue", cex=1.5,
                        main="Pairs plot; correlation between dBhat, dXsq, dDev")
### component + residual plots
        car::crPlots(model,
                     main="Partial (component +) residual plots. Check linear for each predictor.",
                     sub ="xi  vs.  bi*xi + residuals (from full model)")
### Av Plot
        car::avPlots(model, intercept=TRUE,
                     main="Added-variable plots = y~x, adjusted for others, should be linear",
                     cex=1)
### marginal model plots
        par( mfrow=c(nrow1,ncol1), oma = c(0,0,4,0), mar=c(4,6,3,0.5) )
### prevent stopping here due to error
        tryCatch(car::marginalModelPlots(model), error=function(e)e )
        tex1 <- "Marginal model plots (with loess smooth): \n Plots of outcome for each and all predictors"
        graphics::mtext(tex1, line=0, outer = TRUE, cex=1.3)
    }
    par <- p
}
