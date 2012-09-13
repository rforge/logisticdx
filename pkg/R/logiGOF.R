logiGOF <-
function(model){
### require(pROC)
### require(rms)
 source("stukel.R")
 source("logiDx.R")

if (class(model)[1] !="glm") return("LogiGOF only applies to objects of class glm")

### abbreviations:
        f1 <- model
	w1 <- model$data
	nr1 <- nrow(w1)
	m1 <- logiDx(f1)[[1]]

### Pearson residuals by group
PrSj <- sum( m1[m1[,"obs"]>=1,"PeR"]^2) # large S for Sum
DevSj <- sum( m1[m1[,"obs"]>=1,"devR"]^2)

### Pearson residuals by individual subject
PrSi <- sum(residuals(f1, type="pearson")^2)
DevSi <- sum(residuals(f1, type="deviance")^2) # same as Residual Deviance from >summary(f1)

### chi-square tests
### p (by chisq) <0.05 = reject H0, i.e. coefficients are significant predictors

degf1 <- sum(m1[,"obs"]>=1) - (length(f1$coefficients)-1) -1
### degrees freedom = no. covariate patterns (with any observations) - no. predictors in equation +1

### by covariate pattern
ChiPearCov <- list(pValue = 1-pchisq(PrSj, degf1),
	interpret = c("Pearsons Chi-square, calculated by covariate group",
		"p (by chisq) <0.05 = reject H0, i.e. coefficients are significant predictors") )
ChiDevCov <- list(pValue = 1-pchisq(DevSj, degf1),
	interpret = c("Deviance Chi-square, calculated by covariate group",
		"p (by chisq) <0.05 = reject H0, i.e. coefficients are significant predictors") )

### by individual
ChiPearIndiv <- list(pValue = 1-pchisq(PrSi, df.residual(f1)),
	interpret = c("Pearsons Chi-square, calculated by observation",
		"p (by chisq) <0.05 = reject H0, i.e. coefficients are significant predictors") )
ChiDevIndiv <- list(pValue = 1-pchisq(DevSi, df.residual(f1)),
	interpret = c("Deviance Chi-square, calculated by observation",
		"p (by chisq) <0.05 = reject H0, i.e. coefficients are significant predictors") )


### the above doesn't work well as nrow(m1) approaches nrow(w1) so instead use contingency table
m2<- m1[ m1[,"obs"]>=1, c("obs","pi","yhat","y")]
m2 <- cbind(m2,m2[,"obs"]*(1-m2[,"pi"]),m2[,1]-m2[,4] )
colnames(m2) <- c("obs","pi","yhat.y1","y1","yhat.y0","y0") # yNot= probability of y=0 for this pattern
chiFn1 <- function(i,j) sum( ( m2[i,j] - ( m2[i,(j-1)]^2) )/ m2[i,(j-1)] ) # manual chi-sq test
chi1 <- sum(outer(1:nrow(m2), c(4,6), Vectorize(chiFn1)))
degfree1 <- nrow(m2)-ncol(m2)
ChiPearTab <- list(pValue = 1-pchisq(chi1,degfree1),
	interpret=c("Pearsons Chi-square, calculated from table of covariate patterns by outcome",
		"p (by chisq) <0.05 = reject H0, i.e. coefficients are significant predictors") )

GsqFn1 <- function(i,j) 2*sum( m2[i,j] * log(m2[i,j]/m2[i,(j-1)]) )
Gsq1 <- sum( outer(1:nrow(m2), c(4,6), Vectorize(GsqFn1)), na.rm=TRUE )
ChiDevTab <- list(pValue=1-pchisq(Gsq1,degfree1),
	interpret=c("Deviance Chi-square, calculated from table of covariate patterns by outcome",
		"p (by chisq) <0.05 = reject H0, i.e. coefficients are significant predictors") )

### Hosmer Lemeshow GOF test
### p <0.05 reject null hypothesis that the model is a good fit; i.e. model is a poor fit
hosmerlem = function(y, yhat, g=10) {
  cutyhat = cut(yhat,
     breaks = quantile(yhat, probs=seq(0, 1, 1/g)), include.lowest=TRUE)
  obs = xtabs(cbind(1 - y, y) ~ cutyhat)
  expect = xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
  chisq = sum((obs - expect)^2/expect)
  P = 1 - pchisq(chisq, g - 2)
  return(P)
}
### get name of y from formula
HosLem <- list(pValue=hosmerlem(y=w1[,as.character(f1$formula[[2]]) ], yhat=f1$fitted.values, g=10),
	interpret=c("Hosmer & Lemeshow goodness of fit test, with g=10 quantile groups",
		"p <0.05 reject H0 that the model is a good fit; i.e. model is a poor fit",
 		"note may be overly cautions in datasets with large no observations"))

### modified Homser-Lemeshow test
modifiedHL <- function (object, digits = 4) {
    mf <- model.frame(object)
    mf$pred <- fitted(object)
    mf$Residuals <- resid(object)
    mf <- mf[order(mf$pred), ]
    mf$seq <- 1:nrow(mf)
    xcuts <- quantile(1:nrow(mf), prob = seq(0.1:1, by = 0.1))
    mf$group <- cut(mf$seq, xcuts, labels = FALSE)
    mf$group <- ifelse(is.na(mf$group), 0, mf$group)
    mf$group <- as.factor(mf$group)
    ans <- anova(lm(Residuals ~ group, mf))
    pvalue <- lapply(ans, "[[", 1)$"Pr(>F)"
    pval1 <- round(pvalue, digits)
    return(pval1)
}
modHosLem <- list (pValue=modifiedHL(f1),
	interpret=c("modified Hosmer & Lemeshow goodness of fit test, with g=10 quantile groups",
		"p <0.05 reject H0 that the model is a good fit; i.e. model is a poor fit"))

### le Cessie and Houwelingen test
### p<0.05 = reject null hypothesis that the true probabilities are those specified by the model; i.e. model is not a good fit

f2 <- rms::lrm (f1$formula, data=f1$data, x=TRUE, y=TRUE) # All predictors

CesHou <- list(pValue=rms::residuals.lrm(f2,"gof")[[5]],
	interpret=c("le Cessie, van Houwelingen, Copas & Hosmer unweighted sum of squares test for global goodness of fit",
		"p <0.05 reject H0 that the model is a good fit; i.e. model is a poor fit"))

### Osius & Rojek test
### p < 0.05 reject null hypothesis that the true probabilities are those specified by the model; i.e. model is not a good fit
v1 <- m2[,"yhat.y1"]*(1-m2[,"pi"])
c1 <- (1-2*m2[,"pi"] )/ v1
m3<- m1[ m1[,"obs"]>=1, 1:length(f1$coefficients)-1]
lm1 <- lm(c1~ m3, weights=v1) # weighted least squares regression
RSS1 <- sum(lm1$residuals^2)
A1 <- 2*( nrow(m2) - sum(1/m2[,"obs"]) ) # A1 = correction factor for the variance
Xsq1 <- sum(  ( m2[,"y1"] - m2[,"yhat.y1"] )^2 / v1  )
z1 <- ( Xsq1 - ( nrow(m2) - (length(f1$coefficients)-1) -1) )/ sqrt( A1+RSS1)
OsRo <- list (pValue=1-pnorm(z1),
	interpret=c("Osius & Rojek goodness of fit test",
		"significance of Pearson Chi-square by normal approximation (for large samples)",
		"p <0.05 reject H0 that the model is a good fit; i.e. model is a poor fit"))
rm(v1,c1,m3,lm1,RSS1,A1,Xsq1,z1)

### Stukels test
### g1 <- log( f1$fitted.values /(1-f1$fitted.values ) ) # logit for each covariate pattern
### z1 <- 0.5 * g1^2 * ( f1$fitted.values >0.5 )
### z2 <- -0.5 * g1^2 * ( f1$fitted.values <0.5 )
### should use Raos Score test if possible here, otherwise LRT; need recent version of package("stats"), written 2012
### f1a <- update(f1, ~z1)
### stats::anova.glm(f1,f1a, test="LRT")
### f1a <- update(f1, ~z2)
### stats::anova.glm(f1,f1a, test="LRT")

Stuk <- list( pValue=stukel(f1)[3],
	interpret=c("Stukels goodness of fit test",
		"p <0.05 reject H0 that the logistic model is an appropriate link; i.e. consider alternative link"))

### pseudo R^2 tests

### Pearson r^2 correlation of observed outcome with predicted probability
ybar <- sum(m2[,"y1"]) / sum(m2[,"obs"])
pR2a <- (  sum (   (m2[,"y1"] - ( m2[,"obs"] * ybar ) ) * (  m2[,"yhat.y1"] - (m2[,"obs"]*ybar )) ) )^2
pR2b <- (  sum( ( m2[,"y1"] - ybar)^2 )  ) * ( sum(  m2[,"yhat.y1"] - (m2[,"obs"]*ybar )^2 ) )
PR2 <- list("Pearson R^2" = pR2a/pR2b,
	interpret="correlation of observed outcome with predicted")
rm(pR2a,pR2b)

### linear regression-like sum of squares R^2, using covariate patterns
ssR2a <- sum ( (m2[,"y1"]-m2[,"yhat.y1"])^2 )
ssR2b <- sum (  (m2[,"y1"] -( m2[,"obs"]* ybar )  )^2)
ssR2 <- list("sum of squares R^2"=1-(ssR2a/ssR2b), interpret="linear regression-like sum of squares R^2, using covariate patterns")
rm(ssR2a,ssR2b)

### log-likelihood based R^2
fmla <- as.formula(y~1) # intercept-only model
f0 <- glm (fmla, family = binomial("logit"), data=w1) # All predictors
llR2 <- as.numeric(  ( logLik(f0)-logLik(f1) ) /( logLik(f0)- (logLik(f1)+0.5*f1$deviance) )  )
llR2 <- list("Log Likelihood R^2"=llR2,
	interpret="log-likelohood based R^2, calculated by covariate group")

### ROC curve
 e1 <- eval(parse(text="f1$formula[[2]]")) # get name of y/ outcome variable from formula
 e1 <- paste("w1$",e1,sep="")
 e1 <- eval(parse(text=e1)) # change to vector for use in following formula
 roc1 <- pROC::roc (e1 ~ f1$fitted, ci=TRUE, percent=TRUE)[13:14]
 ROC <-list(AUC=roc1[[1]], CI=roc1[[2]])
 rm(e1,roc1)

        res <- list(ChiPearCov=ChiPearCov,
                    ChiPearIndiv=ChiPearIndiv,
                    ChiPearTab=ChiPearTab,
                    OsRo=OsRo,
                    ChiDevCov=ChiDevCov,
                    ChiDevIndiv=ChiDevIndiv,
                    ChiDevTab=ChiDevTab,
                    CovPatTab=m2[,-2],
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
