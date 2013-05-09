##' @name logiDx
##' @export
##' @title Diagnostics for logistic regression
##' @description
##' Returns standard diagnostic measures for a logistic regression model
##'
##' @param model A model of class \bold{glm}
##' @param width if \code{width = TRUE} changes \code{options(width)} to allow
##' all columns in result to be displayed side-by-side
##' @param round1 No. digits to which to round (for display)
##'
##' @return A list containing one matrix, \dQuote{dxMatrix}.
##' \cr \cr
##' The initial columns show all unique combinations of the predictor
##' variables (i.e. one row for each \emph{covariate pattern}).
##' \cr \cr
##' Subsequent columns are labelled as follows:
##'
##' \item{obs}{Number of observations with this covariate pattern}
##' \item{prob}{Probability of this covariate pattern}
##' \item{yhat}{Number of observations of \eqn{y=1}, predicted by the model}
##' \item{y}{\emph{Actual} number of observations of
##' \eqn{y=1} from the data}
##' \item{lev}{\dfn{leverage}, the diagonal of the hat matrix used to
##' generate the model; a measure of influence of this covariate pattern}
##' \item{devR}{\dfn{Deviance residual}, calculated by covariate pattern; a
##' measure of influence of this covariate pattern}
##' \item{PeR}{\dfn{Pearson residual}, calculated by covariate pattern; a
##' measure of influence of this covariate pattern. Given by:
##' \deqn{ \sqrt{obs}\sqrt{\frac{prob}{(1-prob)}}}{
##'  obs^0.5 (prob/1-prob)^0.5}}
##' \item{sPeR}{\dfn{standardized Pearson residual} calculated by covariate
##' pattern; a measure of influence of this covariate pattern. Given by:
##' \deqn{ \frac{PeR}{\sqrt{(1-lev)}}}{
##'  PeR.(1-lev)^0.5}}
##' \item{dBhat}{\dfn{change in Bhat}, the standardized difference between
##' the original maximum likelihood estimates \bold{B} and that the estimates
##' with this covariate pattern excluded}
##' \item{dXsq}{\dfn{change in Chi-square}, decrease in the value of
##' Pearson chi-square statistic with this covariate pattern excluded. Given by:
##' \deqn{sPeR^2}}
##' \item{dDev}{\dfn{change in deviance} \bold{D} with this covariate
##' pattern excluded. Given by:
##' \deqn{ \frac{dev^2}{(1-lev)}}{
##' d^2/(1-lev)}}
##' @note Values for the statistics are calculated by \emph{covariate pattern}.
##' Different values may be obtained if calculated for each individual
##' obervation (i.e. row in data frame).
##' \cr \cr
##' Generally, the values calculated by covariate pattern are preferred,
##' particularly where \eqn{obs >5}.
##' @seealso \code{\link{plotLogiDx}}
##' @keywords array
##' @examples
##' m1 <- genLogiDf()$model
##' logiDx(m1)
##'
logiDx <-
function(model, width=TRUE, round1=3){
    if (class(model)[1] !="glm") return("logiDx only applies to objects of class glm")
    v <- model$coefficients
### no. coefficients, excluding intercept
    l1 <- ifelse(
        any(grepl("Intercept", names(v))),
                 length(v)-1,
                 length(v))
### get data from model
    df1 <- model$data
### make design/ covariate pattern matrix; all combinations of predictors
### creates expression for each unique subset
    lengthUnique <- function(j) paste("length(unique(df1[,",j,"]))",sep="")
    e1t <- sapply(1:l1, lengthUnique)
    e2 <- paste(e1t, collapse="*")
### no. rows = no. covariate patterns
    nr1 <- eval(parse(text=e2))
### no ools; 11x additional measurements
    nc1 <- l1 + 11
### hold results
    m1 <- matrix(NA, ncol=nc1, nrow=nr1)
###
### prob. more efficient method for this section...
###
### make string expression for unique elements in predictors (coefficients)
    Unique <- function(j) paste("unique(df1[,",j,"])", sep="")
### apply this to predictors
    e1t <- sapply(1:l1, Unique)
### convert elements to one string
    e1 <- paste(e1t, collapse=" , ")
### made expression (as string)
    e1 <- paste("expand.grid(",e1,")", sep="")
### convert string to expression
    e1 <- parse(text=e1)
    m1[1:nr1, 1:l1] <- as.matrix(eval(e1))
### arbitrary staring column names
    colnames(m1) <- 1:dim(m1)[2]
### change to coeffient names
    colnames(m1)[1:l1] <- names(v)[2:length(v)]
###
### no. observations per covariate pattern
### string expression:
### subset model data by covariate pattern
    subsetByCov <- function(j) paste0(
        "df1[names(v)[2:length(v)][",j,"]]==m1[i,",j,"]")
    e1t <- sapply(1:l1, subsetByCov)
    e1 <- paste(e1t, collapse=" & ")
### convert string to expression
    e1 <- parse(text=e1)
    subsetByExpr <- function(i) nrow( subset(df1, eval(e1)))
### obs is no observations fitting criteria
    m1[ , l1+1] <- sapply(1:nr1, subsetByExpr)
###
### probability (predicted) each covariate pattern
### creates expression for each unique subset
    makeProb <- function(j) paste0(
        "v[",j+1,"] * m1[i,",j,"]")
    e1t <- sapply(1:l1, makeProb)
    e1 <- paste(e1t, collapse=" + ")
    e1 <- paste("v[1]",e1, sep=" + ")
    e1 <- parse(text=e1)
### eval(e1) gives logit for this covariate pattern
    makeLogit <- function(i) round( exp(eval(e1))/ (1+exp(eval(e1))), 3)
### prob is predicted probability for each pattern
    m1[,l1+2] <- sapply(1:nr1, makeLogit)
###
### yhat =  predicted 'no. outcome=1' in each pattern
### yhat = obs * prob
    m1[ , l1+3] <- m1[ , l1+1] * m1[ , l1+2]
### acutal no.
    e1 <- vector(length = l1, mode="character")
### fun1 <- function(j) paste("df1[xnam[",j,"]]==m1[i,",j,"]", sep="") # creates expression for each unique subsets
# creates expression for each unique subsets
    e1 <- sapply(1:l1, subsetByCov)
    e1 <- paste(e1, collapse=" & ")
    e1 <- parse(text=e1)
### actual no observations fitting these criteria:
    totalY <- function(i) {
	ss <- subset(df1, eval(e1), select="y")
	return( sum(ss$y==1) )
    }
    m1[ ,l1+4] <- sapply(1:nr1, FUN=totalY)
###
    colnames(m1)[(l1+1):(l1+4)] <- c("obs","prob","yhat","y")
###--------------------------------
### make hat matrix
###
### add intercept term to make Design matrix
    d1 <- cbind( rep(1,nrow(m1)), m1[ ,(1:l1)] )
### no. observations fitting criteria
    v <- m1[ , "obs"] * m1[ , "prob"] * (1-(m1[ , "prob"]))
### convert to diagonal matrix and take square root
    v1 <- diag(v)
    v1s <- sqrt(v1)
### hat matrix
    H1 <- v1s %*% d1 %*% solve(t(d1) %*% v1 %*% d1) %*% t(d1) %*% v1s
### hat diagonals = leverage
    m1[ ,l1+5] <- round( diag(H1), round1)
###
### deviance residual for covariate pattern
    devByCov <- function(j){
### y=0 for this covariate pattern
	if (m1[j,"y"] ==0){
            d1 <- log( (1-m1[j, "prob"]) )
            dev <- -sqrt( 2 * m1[j, "obs"] * abs(d1) )
            res <- round(dev, round1)
            return(res)
### y = no. obs. for this covariate pattern
        } else if (m1[j,"y"] == m1[j,"obs"]){
            d1 <- log( (m1[j,"prob"]) )
            dev <- sqrt( 2 * m1[j,"obs"] * abs(d1) )
            res <- round(dev, 3)
            return(res)
        } else {
            d1 <- m1[j,"y"] / (m1[j,"yhat"])
            d2 <- m1[j,"y"] * log(d1)
            d3 <- ( m1[j,"obs"] - m1[j,"y"] ) / ( m1[j,"obs"] * (1 -m1[j,"prob"]) )
            d4 <- (m1[j,"obs"]-m1[j,"y"] )* log(d3)
            d5 <- sqrt( 2*(d2 + d4) )
### 1 if +ve, 0 if -ve
            s1 <- sign(m1[j,"y"]- m1[j,"yhat"])
            dev <- s1*d5
            res <- round(dev, round1)
            return(res)
	}
    }
    m1[ ,l1+6] <- sapply(1:nr1, devByCov)
###
### Pearson residual
    Pear <- function(j){
### y=0 for this covariate pattern
	if (m1[j,"y"] ==0){
            Pr1 <- sqrt( m1[j,"prob"] / (1-m1[j,"prob"]))
            Pr2 <- -sqrt (m1[j,"obs"])
            res <- round( Pr1 * Pr2, 3);res
            return(res)
        } else {
### y>0 for this covariate pattern
            Pr1 <- m1[j,"y"] - m1[j,"yhat"]
            Pr2 <- sqrt(   m1[j,"yhat"] * ( 1-(m1[j,"prob"]) )   )
            res <- round( Pr1/Pr2, 3);res
            return(res)
	}
}
    m1[ ,l1+7] <- sapply(1:nr1, FUN=Pear)
###
    colnames(m1)[(l1+5):(l1+7)] <- c("lev","devR","PeR")
###
### Standardized Pearson residuals
    m1[ ,l1+8] <- round ( m1[ ,"PeR"] / (sqrt (1-m1[ ,"lev"]) ), round1)
    colnames(m1)[(l1+8)] <- "sPeR"
###
### dBhat
### standardized difference in B (maximum likelihood coefficients
### for paramaters) without this pattern
### (standardized by covariance matrix of B);
### should be <1 if little influence on model
###
    m1[ ,l1+9] <- round((m1[ ,"sPeR"])^2*m1[ ,"lev"] / (1-m1[ ,"lev"]), round1)
###
### dXsq
### decrease in Pearson chi-square without this pattern,
### should be <4 if little influence on model
    m1[ ,l1+10] <- round ( (m1[ ,"sPeR"])^2, round1)
###
### dDev
### decrease in deviance without this pattern
    m1[,l1+11] <- round ( (m1[ ,"devR"])^2 / (1-m1[ ,"lev"]) , round1)
###
    colnames(m1)[(l1+9):(l1+11)] <- c("dBhat","dXsq","dDev")
###
    if(isTRUE(width)) options(width= ((l1+11)*6) )
### allows all to be seen together
    res <- list(dxMatrix=m1)
###
### class(res) <- c("logisticDx","matrix")
    return(res)
}
