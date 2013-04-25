logiSS <-
function(model, alpha=0.05, beta=0.8, coeff="x1") {

 if (class(model)[1] !="glm") return("LogiSS only applies to objects of class glm")

 f1 <- model # abbreviate for simplicity

### check if one binomial predictor
 if ( length(f1$coefficients)==2 & length(unique(f1$data[,names(f1$coefficients[2])]))== 2) {
     x1 <- f1$data[,names(f1$coefficients[2])] # convert to 0,1
     x1 <- replace(x1, which(x1== min(x1)), 0)
     x1 <- replace(x1, which(x1== max(x1)), 1)
     P0 <- sum( f1$y[x1==0] ) / length(f1$y)            # P(y=1|x=0)
     P1 <- sum( f1$y[x1==1] ) / length(f1$y)            # P(y=1|x=1)
     Pbar <- P0+P1
     n1 <- (  qnorm(1-alpha)*sqrt(Pbar*(1-Pbar)) + qnorm(beta)*sqrt( P0*(1-P0) + P1*(1-P1) )  )^2 / ((P1-P0)^2)
     tex1 <- paste("Sample size required to show P(y=1|x=0)=",P0," different to P(y=1|x=1)=",P1, sep="")
###alternative method, using sampling distribution of Wald statistic for estimate of logistic regression coefficient
     pi <- sum( x1==0 ) / length(f1$y)
     B1 <- f1$coefficients[[2]]
     W1 <- qnorm(1-alpha)*sqrt( (1/(1-pi)) + (1/pi) )
     W2 <- qnorm(beta)*sqrt( (1/(1-pi)) + (1/(pi*exp(B1))) )
     n2 <- ( (1+2*P0) * (W1+W2)^2 )/ (P0*(B1^2));n2
     tex2 <- paste("Sample size required to show coefficient b1=",B1," rather than b1=0",sep="")
     res <- list(n1=2*n1, interpret1=tex1, n2=n2, interpret2=tex2 )
 }

### check if one continuous predictor
 if ( length(f1$coefficients)==2 & length(unique(f1$data[,names(f1$coefficients[2])]))> 2) {
     x1 <- f1$data[,names(f1$coefficients[2])]
     x1z <- (x1 - mean(x1) )/ sd(x1)    # standardize
     df2 <- data.frame(f1$y,x1z)
     fmla <- as.formula(f1.y~x1z)
     f2 <- glm (fmla, family = binomial("logit"), data=df2) # All predictors
     B0 <- f2$coefficients[[1]]
     P0 <- exp(B0)/( 1+exp(B0) )
     B1 <- f1$coefficients[[2]]
     d1 <- 1+( (1+ (B1^2)) * exp(1.25*(B1^2)) )
     d2 <- 1+exp(-0.25*(B1^2))
     delta <- d1/d2
     W1 <- (  qnorm(1-alpha) + qnorm(beta)*exp(-0.25*(B1^2))  )^2
     n1 <- (1+(2*P0*delta)* (W1/(P0*(B1^2))) )
     res <- list(n1=n1, interpret="Sample size required to show coefficient b1=",B1," rather than b1=0",sep="")
 }

### power size where coefficient is continuous, by method of Hsieh (1989)
 if ( length(f1$coefficients)>2 & length(unique(f1$data[,coeff]))> 2 ) {
     x1 <- f1$data[,coeff]
     x1z <- (x1 - mean(x1) )/ sd(x1)    # standardize
     df2 <- data.frame(f1$y,x1z)
     fmla <- as.formula(f1.y~x1z)
     f2 <- glm (fmla, family = binomial("logit"), data=df2) # All predictors
     B0 <- f2$coefficients[[1]]
     P0 <- exp(B0)/( 1+exp(B0) )
     B1 <- f1$coefficients[names(f1$coefficients)==coeff]
     d1 <- 1+( (1+ (B1^2)) * exp(1.25*(B1^2)) )
     d2 <- 1+exp(-0.25*(B1^2))
     delta <- d1/d2
     W1 <- (  qnorm(1-alpha) + qnorm(beta)*exp(-0.25*(B1^2))  )^2
     n1uc <- (1+(2*P0*delta)* (W1/(P0*(B1^2))) )
### get R^2, multiple correlation between predictor and other predictors in model
     pred1 <- f1$data[ ,!(names(f1$data)==coeff)]
     pred1 <- pred1[,-ncol(pred1)]      # remove y column
     R2 <- summary( lm( f1$y ~ as.matrix(pred1) )) [[8]]
     n1 <- n1uc/(1-R2)
     res <- list (n1=n1,
                  interpret= paste("Sample size required to show coefficient b1=",B1," rather than b1=0",sep=""))
 }

### power size calculation for binary predictor, from Hosmer & Lemeshow
 if ( length(f1$coefficients)>2 & length(unique(f1$data[,coeff])) == 2 ) {

     x1 <- f1$data[,coeff]
### convert to 0,1
     x1 <- replace(x1, which(x1== min(x1)), 0)
     x1 <- replace(x1, which(x1== max(x1)), 1)
     P0 <- sum( f1$y[x1==0] ) / length(f1$y) # P(y=1|x=0)
### method using sampling distribution of Wald statistic for estimate of logistic regression coefficient
### with correction for multiple comparison
     ybar <- sum(f1$y) / length(f1$y)
     PearR2a <- (  sum( (f1$y - ybar)*(f1$fitted.values - ybar) )  )^2
     PearR2b <- sum( (f1$y - ybar)^2) * sum( (f1$fitted.values - ybar)^2 )
     PearR2 <- PearR2a/PearR2b
     rm(PearR2a, PearR2b, ybar)

     pi <- sum( x1==0 ) / length(f1$y)
     B1 <- f1$coefficients[[coeff]]
     W1 <- qnorm(1-alpha)*sqrt( (1/(1-pi)) + (1/pi) )
     W2 <- qnorm(beta)*sqrt( (1/(1-pi)) + (1/(pi*exp(B1))) )
     cor1 <- (1+2*P0)/(1-PearR2)
     n1 <- ( (1+2*P0) * (W1+W2)^2 )/ (P0*(B1^2))
     tex1 <- c( paste("Sample size required to show coefficient b1=",B1," rather than b1=0",sep=""),
               "may be overly conservative; i.e. n may be estimated as too large")
     res <- list(n1=2*n1, interpret=tex1)
 }

### events per covariate
 epc <- min( sum(f1$y), length(f1$y)-sum(f1$y) ) / (length(f1$coefficients)-1)
 epc <- list(epc=epc, interpret="Events per covariate; ideally >10")

### add to result and return
 res1 <- c(res,  epc)
 return(res1)
}
