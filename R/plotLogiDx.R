plotLogiDx <-
function(model, noPerPage=12, identify=FALSE, extras=FALSE,width=1500, height=800) {

    if (class(model)[1] !="glm") return("plotLogiDx only applies to objects of class glm")

### require(pROC)
### require(car)
### require(stats)
### require(graphics)

### abbreviations
    m1 <- logiDx(model)$dxMatrix
    f1 <- model
    w1 <- model$data

fun1 <- function(x) (x^2 + noPerPage )/ x
nrow1 <- round( optimize(fun1, interval=seq(1:noPerPage))$minimum, 0)
ncol1 <- round( noPerPage/nrow1, 0); nrow1;ncol1

windows(record=T, width=width, height=height)
p <- par

par(mfrow= c(nrow1,ncol1),oma = c(0, 0, 4, 0), mar=c(4,6,3,0.5))
###oma = outer margins, mar=margins, bottom,left,top,right

### colours()[grep("blue",colours())] # select colours
###-------------------------------------------
### probabiility by leverage

graphics::plot(m1[,"pi"], m1[,"lev"], xlab="Probability for this covariate pattern", main="Probability by leverage", ylab="Leverage (hat matrix diagonal)")
###graphics::mtext("Leverage \n (hat \n matrix \n diagonal)", side=2, line=3,las=1) # horizontal on y-axis
points(m1[,"pi"], m1[,"lev"], pch=21, cex = 2, col = c("deepskyblue","dodgerblue"), bg=c("deepskyblue","dodgerblue"))
if (isTRUE(identify)) graphics::identify(m1[,"pi"], m1[,"lev"])

txt <- paste("Diagnostic plots for logistic regression \n", deparse(f1$formula), sep=" ")
graphics::mtext(txt, line = 0.3, outer = TRUE)

###-------------------------------------------
### probability by ...

graphics::plot(m1[,"pi"], m1[,"dXsq"], xlab="Probability for this covariate pattern", ylab="dXsq = decrease in Pearson Chi-sq \n without this pattern", main="dXsq by probability")
points(m1[,"pi"], m1[,"dXsq"], pch=21, cex = 2, col = c("deepskyblue","dodgerblue"), bg=c("deepskyblue","dodgerblue"))
if (isTRUE(identify)) graphics::identify(m1[,"pi"], m1[,"dXsq"])

txt<- paste("Diagnostic plots for logistic regression \n", deparse(f1$formula), sep=" ")
graphics::mtext(txt, line = 0.3, outer = TRUE)

graphics::plot(m1[,"pi"], m1[,"dBhat"], xlab="Probability for this covariate pattern", ylab="dBhat = decrease in Bhat \n without this pattern", main="dBhat by probability")
points(m1[,"pi"], m1[,"dBhat"], pch=21, cex = 2, col = c("deepskyblue","dodgerblue"), bg=c("deepskyblue","dodgerblue"))
if (isTRUE(identify)) graphics::identify(m1[,"pi"], m1[,"dBhat"])

txt <- paste("Diagnostic plots for logistic regression \n", deparse(f1$formula), sep=" ")
graphics::mtext(txt, line = 0.3, outer = TRUE)

graphics::plot(m1[,"pi"], m1[,"dDev"], xlab="Probability for this covariate pattern", ylab="dDev = decrease in Deviance \n without this pattern", main="dDev by probability")
points(m1[,"pi"], m1[,"dDev"], pch=21, cex = 2, col = c("deepskyblue","dodgerblue"), bg=c("deepskyblue","dodgerblue"))
if (isTRUE(identify)) graphics::identify(m1[,"pi"], m1[,"dDev"])

###---------------------------
### bubble plot

radius <- sqrt(m1[,"dBhat"]/ m1[,"pi"])
symbols(m1[,"pi"], m1[,"dXsq"], circles=radius, inches=0.35, fg="white", bg=c("deepskyblue","dodgerblue"), xlab="Probability for this covariate pattern", ylab="dXsq = decrease in Pearson Chi-sq \n without this pattern", main = "Area proportional to dBhat \n Decrease in Bhat without this pattern")
if (isTRUE(identify)) graphics::identify(m1[,"pi"], m1[,"dXsq"])

txt<- paste("Diagnostic plots for logistic regression \n", deparse(f1$formula), sep=" ")
graphics::mtext(txt, line = 0.3, outer = TRUE)

###-------------------------------------------
### leverage by...

graphics::plot(m1[,"lev"], m1[,"dXsq"], xlab="Leverage for this covariate pattern", ylab="dXsq = decrease in Pearson Chi-sq \n without this pattern", main="dXsq by leverage")
points(m1[,"lev"], m1[,"dXsq"], pch=21, cex = 2, col = c("deepskyblue","dodgerblue"), bg=c("deepskyblue","dodgerblue"))
if (isTRUE(identify)) graphics::identify(m1[,"lev"], m1[,"dXsq"])

txt <- paste("Diagnostic plots for logistic regression \n", deparse(f1$formula), sep=" ")
graphics::mtext(txt, line = 0.3, outer = TRUE)

graphics::plot(m1[,"lev"], m1[,"dBhat"], xlab="Leverage for this covariate pattern", ylab="dBhat = decrease in Bhat \n without this pattern", main="dBhat by leverage")
points(m1[,"lev"], m1[,"dBhat"], pch=21, cex = 2, col = c("deepskyblue","dodgerblue"), bg=c("deepskyblue","dodgerblue"))
if (isTRUE(identify)) graphics::identify(m1[,"lev"], m1[,"dBhat"])

txt<- paste("Diagnostic plots for logistic regression \n", deparse(f1$formula), sep=" ")
graphics::mtext(txt, line = 0.3, outer = TRUE)

graphics::plot(m1[,"lev"], m1[,"dDev"], xlab="Leverage for this covariate pattern", ylab="dDev = decrease in Deviance \n without this pattern", main="dDev by leverage")
points(m1[,"lev"], m1[,"dDev"], pch=21, cex = 2, col = c("deepskyblue","dodgerblue"), bg=c("deepskyblue","dodgerblue"))
if (isTRUE(identify)) graphics::identify(m1[,"lev"], m1[,"dDev"])

txt <- paste("Diagnostic plots for logistic regression \n", deparse(f1$formula), sep=" ")
graphics::mtext(txt, line = 0.3, outer = TRUE)

### ROC curve
e1 <- eval(parse(text="f1$formula[[2]]")) # get name of y/ outcome variable from formula
e1 <- paste("w1$",e1,sep="")
e1 <- eval(parse(text=e1)) # change to vector for use in following formula
r1 <- pROC::roc (e1 ~ f1$fitted, ci=TRUE, percent=TRUE)
pROC::plot.roc(r1, print.auc=TRUE, grid=TRUE, print.auc.cex=0.8, main="ROC curve") # 0.5 = chance, aim >0.7

txt <- paste("Diagnostic plots for logistic regression \n", deparse(f1$formula), sep=" ")
graphics::mtext(txt, line = 0.3, outer = TRUE)

### influence plot
if (isTRUE(identify)){
car::influencePlot(f1, id.n=1.5, id.method="identify",id.col="blue", scale=15,
	xlab="Hat values (vertical lines at x2, x3 average hat value)",
	main="Area proportional to Cooks distance",
	sub="Cooks = change in coefficients if this point dropped")
	}else{
	car::influencePlot(f1, id.n=1.5, id.method="noteworthy",id.col="blue", scale=15,
	xlab="Hat values (vertical lines at x2, x3 average hat value)",
	main="Area proportional to Cooks distance",
	sub="Cooks = change in coefficients if this point dropped")
	}

txt <- paste("Diagnostic plots for logistic regression \n", deparse(f1$formula), sep=" ")
graphics::mtext(txt, line = 0.3, outer = TRUE)


###---------------------------------------------------------
### optional plots

if (isTRUE(extras)) {
### studentized residuals vs hat values
    graphics::plot(rstudent(f1) ~ hatvalues(f1),
         xlab="Hat values",
         ylab="Studentized residuals (residual / estimate of its S.D.)",
         main="Studentized residuals by hat values") # recommended by some
    points(hatvalues(f1), rstudent(f1), pch=21, cex = 2, col = c("deepskyblue","dodgerblue"), bg= c("deepskyblue","dodgerblue"))
    if (isTRUE(identify)) {graphics::identify(hatvalues(f1), rstudent(f1), cex=1.5)}
    txt <- paste("Diagnostic plots for logistic regression \n", deparse(f1$formula), sep=" ")
    graphics::mtext(txt, line = 0.3, outer = TRUE)

### spread level plot
   car::spreadLevelPlot(f1, layout=NA,
                    main="Spread-Level plot \n Spread should be even \n if variance of error is constant")
    txt <- paste("Diagnostic plots for logistic regression \n", deparse(f1$formula), sep=" ")
    graphics::mtext(txt, line = 0.3, outer = TRUE)

### qqnorm for residuals
    q1 <- stats::qqnorm(statmod::qresid(f1),col="darkblue", cex=1.5,
                 main="Normal Q-Q plot for residuals from formula \n Check if residuals normally distributed",
                 ylab="Quantiles from formula")
    stats::qqline((statmod::qresid(f1)),col="darkblue",lwd=1.5)
    if (isTRUE(identify)) graphics::identify(q1, cex=1.5)
    txt <- paste("Diagnostic plots for logistic regression \n", deparse(f1$formula), sep=" ")
    graphics::mtext(txt, line = 0.3, outer = TRUE)

### Cooks, stud~ res~ and hat values by index
    if (isTRUE(identify)){
	car::influenceIndexPlot(f1, col="dodgerblue", id.method = "identify",
                           vars=c("Cook", "Studentized", "hat"),
                           main="Cooks distance, studentized residuals \n and hat values by Index (of observations)")
    } else {
        car::influenceIndexPlot(f1, col="dodgerblue", id.method = "y",
                           vars=c("Cook", "Studentized", "hat"),
                           main="Cooks distance, studentized residuals \n and hat values by Index (of observations)")
    }

### pairs plot for measures of influence
    graphics::pairs(m1[,c("dBhat", "dXsq", "dDev")], col="darkblue", cex=1.5,
          main="Pairs plot; correlation between dBhat, dXsq, dDev")

### component + residual plots
    car::crPlots(f1,
            main="Partial (component +) residual plots. Check linear for each predictor.",
            sub ="xi  vs.  bi*xi + residuals (from full model)")

### Av Plot
    car::avPlots(f1, intercept=TRUE,
            main="Added-variable plots = y~x, adjusted for others, should be linear", cex=1)

### marginal model plots
par(mfrow= c(nrow1,ncol1),oma = c(0, 0, 4, 0), mar=c(4,6,3,0.5))
    car::marginalModelPlots(model=f1)
txt <- "Marginal model plots: \n Plots of outcome for each and all predictors \n (with loess smooth)"
graphics::mtext(txt, line=0, outer = FALSE, cex=1.3)


}

par <- p

}
