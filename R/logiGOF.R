##' @name logiGOF
##' @export
##' @include stukel.R
##' @include logiDx.R
##' @title Goodness of fit tests for a logistic regression model
##' @description
##' Gives 15 commonly employed measures of goodness of fit for a logistic
##' regression model
##'
##' @param model A model of class \bold{glm}
##' @param g No. groups (quantiles) into which to split observatiosn for Hosmer-Lemeshow and modified Hosmer-Lemeshow tests.
##' @return A list with the following items:
##' \item{ChiPearCov}{Pearsons Chi-square, calculated by covariate group,
##' with p value and interpretation}
##' \item{ChiPearIndiv}{Pearsons Chi-square, calculated by individual
##' observation, with p value and interpretation}
##' \item{ChiPearTab}{Pearsons Chi-square, calculated by table of covariate
##' patterns by outcome, with p value and interpretation}
##' \item{OsRo}{Osius & Rojek test of the logistic link, with p value and
##' interpretation}
##' \item{ChiDevCov}{Deviance Chi-square, calculated by covariate group,
##' with p value and interpretation}
##' \item{ChiDevIndiv}{Deviance Chi-square, calculated by individual
##' observation, with p value and interpretation}
##' \item{ChiDevTab}{Deviance Chi-square, calculated by table of covariate
##' patterns by outcome, with p value and interpretation}
##' \item{CovPatTab}{Matrix of covariance patterns, used to calculate above
##' Chi-square tests of Pearson residuals and Deviance}
##' \item{HosLem}{Hosmer & Lemeshow goodness of fit test, with \emph{g=10}
##' quantile groups,with p value and interpretation}
##' \item{modHosLem}{modified Hosmer & Lemeshow goodness of fit test, with
##' \emph{g=10} quantile groups, with p value and interpretation}
##' \item{CesHou}{le Cessie, van Houwelingen, Copas & Hosmer unweighted sum
##' of squares test for global goodness of fit, with p value and interpretation}
##' \item{Stuk}{Stukels test of the appropriateness of the logistic link,
##' with p value and interpretation}
##' \item{PR2}{Pearsons R^2, correlation of observed outcome with predicted}
##' \item{ssR2}{linear regression-like sum of squares R^2, using covariate
##' patterns}
##' \item{llR2}{log-likelohood based R^2, calculated by covariate group}
##' \item{ROC}{Area under the Receiver Operating Curve, with 95\% CI by
##' method of DeLong}
##' @note Warning: Will fail if cannot generate a hat matrix for the model
##' using \code{logiDx}
##'
##' @author Yongmei Ni: modified Hosmer & Lemeshow goodness of fit test (adapted)
##' @seealso \code{\link{logiDx}}
##' @keywords htest
##' @examples
##'
##' set.seed(1)
##' m1 <- genLogiDf(n=100)$model
##'
logiGOF <- function(model, g=10){
### check if model is glm
    err1 <- "LogiGOF only applies to objects of class glm"
    if (class(model)[1] !="glm") stop(err1)
### abbreviations:
    df1 <- model$data
    m1 <- logiDx(model)[[1]]
### Pearson residuals by group
### (large S for Sum)
    PrSj <- sum( m1[ m1[ ,"obs"]>=1, "PeR" ]^2 )
    DevSj <- sum( m1[ m1[ ,"obs"]>=1, "devR" ]^2 )
###
### Pearson residuals by individual subject
    PrSi <- sum(residuals(model, type="pearson")^2)
### same as Residual Deviance from >summary(model)
    DevSi <- sum(residuals(model, type="deviance")^2)
### chi-square tests
### p (by chisq) <0.05 = reject H0,
### i.e. coefficients are significant predictors
###
### degrees freedom = no. covariate patterns (with any observations) - no. predictors in equation +1
    degf1 <- sum(m1[ ,"obs"]>=1) - (length(model$coefficients)-1) -1
### by covariate pattern
    ChiPearCov <- list(pValue = 1-pchisq(PrSj, degf1),
                       interpret = c("Pearsons Chi-square, calculated by covariate group",
                       "p (by chisq) <0.05 = reject H0, i.e. coefficients are significant predictors") )
    ChiDevCov <- list(pValue = 1-pchisq(DevSj, degf1),
                      interpret = c("Deviance Chi-square, calculated by covariate group",
                      "p (by chisq) <0.05 = reject H0, i.e. coefficients are significant predictors") )
###
### by individual
    ChiPearIndiv <- list(pValue = 1-pchisq(PrSi, df.residual(model)),
                         interpret = c("Pearsons Chi-square, calculated by observation",
                         "p (by chisq) <0.05 = reject H0, i.e. coefficients are significant predictors") )
    ChiDevIndiv <- list(pValue = 1-pchisq(DevSi, df.residual(model)),
                        interpret = c("Deviance Chi-square, calculated by observation",
                        "p (by chisq) <0.05 = reject H0, i.e. coefficients are significant predictors") )
###
### the above doesn't work well as nrow(m1) approaches nrow(model$data) so instead use contingency table
    m2<- m1[ m1[ ,"obs"]>=1, c("obs","prob","yhat","y")]
    m2 <- cbind(m2, m2[ ,"obs"] * ( 1-m2[ ,"prob"]), m2[ ,1] -m2[ ,4] )
    colnames(m2) <- c("obs","prob","yhat.y1","y1","yhat.y0","y0")
###
### manual chi-sq test
    chiFn1 <- function(i,j) sum( ( m2[i,j] - ( m2[i,(j-1)]^2) )/ m2[i,(j-1)] )
    chi1 <- sum( outer(1:nrow(m2), c(4,6), Vectorize(chiFn1)) )
    degfree1 <- nrow(m2)-ncol(m2)
    ChiPearTab <- list(pValue = 1-pchisq(chi1,degfree1),
                       interpret=c("Pearsons Chi-square, calculated from table of covariate patterns by outcome",
                       "p (by chisq) <0.05 = reject H0, i.e. coefficients are significant predictors") )
###
    GsqFn1 <- function(i,j) 2*sum( m2[i,j] * log(m2[i,j]/m2[i,(j-1)]) )
    Gsq1 <- sum( outer(1:nrow(m2), c(4,6), Vectorize(GsqFn1)), na.rm=TRUE )
    ChiDevTab <- list(pValue = 1-pchisq(Gsq1,degfree1),
                      interpret =c("Deviance Chi-square, calculated from table of covariate patterns by outcome",
                      "p (by chisq) <0.05 = reject H0, i.e. coefficients are significant predictors") )
###
### Hosmer Lemeshow GOF test
### p <0.05 reject null hypothesis that the model is a good fit; i.e. model is a poor fit
    hosmerLem <- function(y, yhat, g1=g) {
        cutyhat <- cut(yhat,
                       breaks = quantile(yhat, probs=seq(0, 1, 1/g1)),
                       include.lowest=TRUE)
        obs <- xtabs(cbind(1 - y, y) ~ cutyhat)
        expect <- xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
        chisq <- sum((obs - expect)^2/expect)
        P <- 1 - pchisq(chisq, g - 2)
        return(P)
    }
###
    HosLem <- list(pValue=hosmerLem(
                   y = model$data[[ncol(model$data)]],
                   yhat = model$fitted.values,
                   g = 10),
                   interpret=c(
                   paste0("modified Hosmer & Lemeshow goodness of fit test, with g=",g," quantile groups"),
                   "p <0.05 reject H0 that the model is a good fit; i.e. model is a poor fit",
                   "note may be overly cautions in datasets with large no observations"))
###
### modified Homser-Lemeshow test
    modHL <- function (model, g1=g) {
        xcuts <- stats::quantile(1:nrow(model$data), prob = seq(0.1:1, by=1/g1))
        group <- cut(1:nrow(model$data), xcuts, labels = FALSE)
        group <- ifelse(is.na(group), 0, group)
        group <- as.factor(group)
        Residuals <- resid(model)[order(fitted(model))]
        ans <- stats::anova(lm(Residuals ~ group))
        pvalue <- lapply(ans, "[[", 1)$"Pr(>F)"
        return(pvalue)
    }
    modHosLem <- list (pValue=modHL(model),
                       interpret=c(
                       paste0("modified Hosmer & Lemeshow goodness of fit test, with g=",g," quantile groups"),
                       "p <0.05 reject H0 that the model is a good fit; i.e. model is a poor fit"))
###
### le Cessie and Houwelingen test
### p<0.05 = reject null hypothesis that the true probabilities are those specified by the model; i.e. model is not a good fit
### All predictors
    f2 <- rms::lrm (model$formula, data=model$data, x=TRUE, y=TRUE)
    CesHou <- list(pValue=rms::residuals.lrm(f2,"gof")[[5]],
                   interpret=c("le Cessie, van Houwelingen, Copas & Hosmer unweighted sum of squares test for global goodness of fit",
                   "p <0.05 reject H0 that the model is a good fit; i.e. model is a poor fit"))
###
### Osius & Rojek test
### p < 0.05 reject null hypothesis that the true probabilities are those specified by the model; i.e. model is not a good fit
    v1 <- m2[ ,"yhat.y1"]*(1-m2[ ,"prob"])
    c1 <- (1-2*m2[ ,"prob"] )/ v1
    m3<- m1[ m1[ ,"obs"]>=1, 1:length(model$coefficients)-1]
### weighted least squares regression
    lm1 <- stats::lm(c1~ m3, weights=v1)
    RSS1 <- sum(lm1$residuals^2)
### A1 = correction factor for the variance
    A1 <- 2*( nrow(m2) - sum(1/m2[,"obs"]) )
    Xsq1 <- sum(  ( m2[ ,"y1"] - m2[ ,"yhat.y1"] )^2 / v1  )
    z1 <- ( Xsq1 - ( nrow(m2) - (length(model$coefficients)-1) -1) )/ sqrt( A1+RSS1)
    OsRo <- list (pValue=1-pnorm(z1),
                  interpret=c("Osius & Rojek goodness of fit test",
                  "significance of Pearson Chi-square by normal approximation (for large samples)",
                  "p <0.05 reject H0 that the model is a good fit; i.e. model is a poor fit"))
### Stukels test
###
### g1 <- log( model$fitted.values /(1-model$fitted.values ) ) # logit for each covariate pattern
### z1 <- 0.5 * g1^2 * ( model$fitted.values >0.5 )
### z2 <- -0.5 * g1^2 * ( model$fitted.values <0.5 )
### should use Raos Score test if possible here, otherwise LRT; need recent version of package("stats"), written 2012
### f1a <- update(model, ~z1)
### stats::anova.glm(model,f1a, test="LRT")
### f1a <- update(model, ~z2)
### stats::anova.glm(model,f1a, test="LRT")
    Stuk <- list( pValue=stukel(model)[3],
                 interpret=c("Stukels goodness of fit test",
                 "p <0.05 reject H0 that the logistic model is an appropriate link; i.e. consider alternative link"))
###
### pseudo R^2 tests
###
### Pearson r^2 correlation of observed outcome with predicted probability
    ybar <- sum(m2[ ,"y1"]) / sum(m2[ ,"obs"])
    pR2a <- (  sum(   (m2[ ,"y1"] - ( m2[ ,"obs"] * ybar ) ) *
                   (  m2[ ,"yhat.y1"] - (m2[ ,"obs"]*ybar )) )  )^2
    pR2b <- (  sum( ( m2[ ,"y1"] - ybar)^2 )  ) *
        ( sum(  m2[ ,"yhat.y1"] - (m2[ ,"obs"]*ybar )^2 ) )
    PR2 <- list("Pearson R^2" = pR2a/pR2b,
                interpret="correlation of observed outcome with predicted")
###
### linear regression-like sum of squares R^2, using covariate patterns
    ssR2a <- sum ( (m2[ ,"y1"]-m2[ ,"yhat.y1"])^2 )
    ssR2b <- sum (  (m2[ ,"y1"] -( m2[ ,"obs"]* ybar )  )^2)
    ssR2 <- list("sum of squares R^2"=1-(ssR2a/ssR2b),
                 interpret="linear regression-like sum of squares R^2, using covariate patterns")
###
### log-likelihood based R^2
###
### intercept-only model
    fmla <- as.formula(y~1)
    f0 <- glm (fmla, family = binomial("logit"), data=model$data)
### use index [1] to return numeric value
    llR2 <- ( ( stats::logLik(f0) -stats::logLik(model) ) /
             (stats::logLik(f0)- (stats::logLik(model)+0.5*model$deviance) ))[1]
    llR2 <- list("Log Likelihood R^2"=llR2,
                 interpret="log-likelohood based R^2, calculated by covariate group")
###
### ROC curve
###
### y/outcome variable from formula is last column in model$data
### [13:14] gets area under curve and confidence interval from roc
    roc1 <- pROC::roc (model$data[ncol(model$data)][[1]] ~ model$fitted,
                       ci=TRUE, percent=TRUE)[13:14]
    ROC <-list(AUC=roc1[[1]],
               CI=roc1[[2]])
###
    res <- list(ChiPearCov=ChiPearCov,
                ChiPearIndiv=ChiPearIndiv,
                ChiPearTab=ChiPearTab,
                OsRo=OsRo,
                ChiDevCov=ChiDevCov,
                ChiDevIndiv=ChiDevIndiv,
                ChiDevTab=ChiDevTab,
                CovPatTab=m2[ ,-2],
                HosLem=HosLem,
                modHosLem=modHosLem,
                CesHou=CesHou,
                Stuk=Stuk,
                PR2=PR2,
                ssR2=ssR2,
                llR2=llR2,
                ROC=ROC)
    return(res)
}
