##' @name logiSS
##' @export
##' @title Sample size for given coefficient and events per covariate for model
##'
##' @description
##' Gives sample size necessary to demonstrate that coefficient in model for
##' given predictor is equal to its given value (rather than equal to zero)
##' for a given level of power and significance.
##' \cr \cr
##' Also number of events (smaller of outcome \eqn{y=0}
##' and outcome \eqn{y=1}) per predictor.
##' \cr \cr
##' Uses different methods depending on whether model has one binomial, one
##' continuous or multiple predictors.
##'
##' @param model A logistic regression model of class \bold{glm}
##' @param alpha significance level \eqn{\alpha}
##' for null-hypothesis significance test
##' @param beta power \eqn{\beta} for null-hypothesis significance test
##' @param coeff Name of predictor (coefficient) in model to be tested
##' @return A list of:
##' \item{res}{Result: Sample size required to show coefficient for
##' predictor is as given in the model rather than 0}
##' \item{epc}{Events per covariate; should be >10 to make meaningful
##' statements about coefficients obtained}
##' @keywords htest
##' @examples
##' set.seed(1)
##' ### one coefficient, which is binomial
##' f1 <- genLogiDf(b=1,c=0,n=50)$model
##' logiSS(f1)
##' ###
##' ### one coefficient, which is continuous
##' f1 <- genLogiDf(f=0,b=0,c=1,n=50,asFactor=TRUE)$model
##' logiSS(f1, coeff="x1")
##' ###
##' ### binomail coefficient
##' f1 <- genLogiDf(f=1,b=1,c=1,n=50)$model
##' logiSS(f1,coeff="x2")
##' ###
##' ### continuous coefficient
##' f1 <- genLogiDf(f=1,b=1,c=1,n=50)$model
##' logiSS(f1,coeff="x3")
##'
logiSS <- function(model, alpha=0.05, beta=0.8, coeff="x1") {
    if (class(model)[1] !="glm")
        stop("LogiSS only applies to objects of class glm")
    if (! coeff %in% names(model$coefficients))
        stop(paste0("Coefficient ",coeff," must be in model"))
### check if one binomial coefficient
 if ( length(model$coefficients)==2 &
     length(unique(model$data[,names(model$coefficients[2])]))== 2 ) {
### convert to 0,1
     x1 <- model$data[,names(model$coefficients[2])]
     x1 <- replace(x1, which(x1== min(x1)), 0)
     x1 <- replace(x1, which(x1== max(x1)), 1)
### probability P(y=1|x=0)
     P0 <- sum( model$y[x1==0] ) / length(model$y)
### P(y=1|x=1)
     P1 <- sum( model$y[x1==1] ) / length(model$y)
     Pbar <- P0+P1
     n1 <- (  qnorm(1-alpha)*sqrt(Pbar*(1-Pbar)) +
            qnorm(beta)*sqrt( P0*(1-P0) + P1*(1-P1) )  )^2 / ((P1-P0)^2)
     tex1 <- paste0("Sample size required to show P(y=1|x=0) =",
                    P0, " different to P(y=1|x=1) =", P1)
### alternative method, using sampling distribution of Wald statistic
### for estimate of logistic regression coefficient
     pi <- sum( x1==0 ) / length(model$y)
### coefficient for predictor
     B1 <- model$coefficients[[2]]
     W1 <- stats::qnorm(1-alpha)*sqrt( (1/(1-pi)) + (1/pi) )
     W2 <- stats::qnorm(beta)*sqrt( (1/(1-pi)) + (1/(pi*exp(B1))) )
     n2 <- ( (1+2*P0) * (W1+W2)^2 )/ (P0*(B1^2))
     tex2 <- paste0("Sample size required to show coefficient b1 =",
                    B1, " rather than b1 =0")
     res <- list(n1=2*n1, interpret1=tex1, n2=n2, interpret2=tex2 )
 }
###
### check if one continuous coefficient (>2 unique values)
    if ( length(model$coefficients)==2 &
        length(unique(model$data[,names(model$coefficients[2])]))> 2) {
        x1 <- model$data[,names(model$coefficients[2])]
### standardize
        x1z <- (x1 - mean(x1) )/ sd(x1)
        f2 <- glm (model$y ~ x1z, family = binomial("logit"))
        B0 <- f2$coefficients[[1]]
### convert to probability
        P0 <- exp(B0)/( 1+exp(B0) )
        B1 <- model$coefficients[[2]]
        d1 <- 1+( (1+ (B1^2)) * exp(1.25*(B1^2)) )
        d2 <- 1+exp(-0.25*(B1^2))
        delta <- d1/d2
        W1 <- (  stats::qnorm(1-alpha) +
               stats::qnorm(beta)*exp(-0.25*(B1^2))  )^2
        n1 <- (1+(2*P0*delta)* (W1/(P0*(B1^2))) )
        tex1 <- paste0("Sample size required to show coefficient b1= ", B1," rather than b1 =0",sep="")
        res <- list(n1=n1, interpret=tex1)
     }
###
### more than one coefficient:
###
### power size where coefficient is continuous, by method of Hsieh (1989)
    if ( length(model$coefficients)>2 &
        length(unique(model$data[,coeff]))> 2){
        x1 <- model$data[,coeff]
        x1z <- (x1 - mean(x1) )/ sd(x1)
        f2 <- glm (model$y ~ x1z, family = binomial("logit"))
        B0 <- f2$coefficients[[1]]
        P0 <- exp(B0)/( 1+exp(B0) )
        B1 <- model$coefficients[names(model$coefficients)==coeff]
        d1 <- 1+( (1+ (B1^2)) * exp(1.25*(B1^2)) )
        d2 <- 1+exp(-0.25*(B1^2))
        delta <- d1/d2
        W1 <- (  stats::qnorm(1-alpha) + stats::qnorm(beta)*exp(-0.25*(B1^2))  )^2
        n1uc <- (1+(2*P0*delta)* (W1/(P0*(B1^2))) )
### get R^2, multiple correlation between predictor and other predictors in model
        pred1 <- model$data[ ,!(names(model$data)==coeff)]
### remove y column
        pred1 <- pred1[,-ncol(pred1)]
        R2 <- summary( lm( model$y ~ as.matrix(pred1) )) [[8]]
        n1 <- n1uc/(1-R2)
        tex1 <- paste0("Sample size required to show coefficient ",
                       coeff," =",B1," rather than ",coeff," =0")
        res <- list (n1=n1, interpret= tex1)
 }
###
### power size calculation for binary predictor, from Hosmer & Lemeshow
    if ( length(model$coefficients)>2 &
        length(unique(model$data[,coeff])) == 2 ){
        x1 <- model$data[,coeff]
### convert to 0,1
        x1 <- replace(x1, which(x1== min(x1)), 0)
        x1 <- replace(x1, which(x1== max(x1)), 1)
### P(y=1|x=0)
        P0 <- sum( model$y[x1==0] ) / length(model$y)
### method using sampling distribution of Wald statistic
### for estimate of logistic regression coefficient
### with correction for multiple comparison
        ybar <- sum(model$y) / length(model$y)
### numerator
        PearR2n <- (  sum( (model$y - ybar)*(model$fitted.values - ybar) )  )^2
### denominator
        PearR2d <- sum( (model$y - ybar)^2) * sum( (model$fitted.values - ybar)^2 )
        PearR2 <- PearR2n/PearR2d
        pi <- sum( x1==0 ) / length(model$y)
        B1 <- model$coefficients[names(model$coefficients)==coeff]
        W1 <- stats::qnorm(1-alpha)*sqrt( (1/(1-pi)) + (1/pi) )
        W2 <- stats::qnorm(beta)*sqrt( (1/(1-pi)) + (1/(pi*exp(B1))) )
        cor1 <- (1+2*P0)/(1-PearR2)
        n1 <- ( (1+2*P0) * (W1+W2)^2 )/ (P0*(B1^2))
        tex1 <- c(paste0("Sample size required to show coefficient ",
                       coeff," =",B1," rather than ",coeff," =0"),
                  "may be overly conservative; i.e. n may be estimated as too large")
     res <- list(n1=2*n1, interpret=tex1)
 }
### events per covariate
    epc <- min( sum(model$y), length(model$y)-sum(model$y) ) / (length(model$coefficients)-1)
    epc <- list(epc=epc, interpret="Events per covariate; ideally >10")
###
### add to result and return
    res1 <- c(res,  epc)
    return(res1)
}
