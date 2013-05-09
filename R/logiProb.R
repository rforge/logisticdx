##' @name logiProb
##' @export
##' @title Logits, odds ratios and probabilities for all
##' combinations of coefficients
##' @description
##' For all combinations of predictors in model (including intercept only),
##' generate logit, log odds ratio (OR), OR and probability.
##'
##' Values are calculated for a change in the value of the coeffient for the
##' predictor from 0 to 1. (For continuous predictors changes of more than one
##' unit may have more practical significance).
##'
##' @param model A logistic regression model of class \bold{glm}
##' @return A matrix giving, for each combination of predictors (coefficients):
##' \item{logit}{The logit for a given combination of coefficients}
##' \item{ln(OR)}{Natural log of odds ratio}
##' \item{OR}{Odds ratio}
##' \item{p}{probability}
##' @keywords array math
##' @examples
##' set.seed(1)
##' f1 <- genLogiDf(n=50)$model
##' logiProb(f1)
##'
logiProb <- function(model){
### require(gRbase)
### for gRbase::combnPrim instead of utils::combn (works faster)
### the following dependencies may be necessary, install as follows:
### source("http://bioconductor.org/biocLite.R")
### biocLite("graph")
### biocLite("RBGL")
    if (class(model)[1] !="glm") stop("LogiProb only applies to objects of class glm")
###
    v <- model$coefficients
### short for intercept
    names(v)[1] <- "Int"
### remove intercept
    v1 <- v[-1]
### used to generate combinations without intercept
### binomial expansion; gives no. of logits
    ch1 <- choose(n=length(v1), k=0:length(v1))
### hold results
### no. rows = no. combinations of coefficients
    res <- matrix(NA, nrow=sum(ch1), ncol=4,
                  dimnames = list(1:sum(ch1),c("logit","ln(OR)","OR","p")))
### remove first and last values from ch1:
    ch2 <- ch1[-c(1,length(ch1))]
### names
### generate combinations of coefficient names
    genComCoefName <- function(i){
### all combinations of coefficients with 'i' elments
        co1 <- gRbase::combnPrim(names(v1),i)
        combnPrim(names(v1),i)
### change to character vector form
        pasteAddCol <- function(j) paste(co1[,j],collapse="+")
        co1 <- sapply(1:ncol(co1), FUN=pasteAddCol, USE.NAMES=FALSE)
### add intercept term to character vector
        addInt <- function(k) paste (names(v)[1], co1[k], sep="+")
        co1 <- sapply(1:length(co1), FUN=addInt)
        return(co1)
    }
### first element (intercept only)
    rownames(res)[1] <- "Int"
### middle elemnts
    rownames(res)[-c(1,nrow(res))]  <- unlist(
        sapply(1:(length(ch2)),FUN=genComCoefName))
### last element (intercept + all coefficients)
    rownames(res)[nrow(res)] <- paste(names(v), collapse="+")
###
### generate logits (from combinations of coefficients)
    genLogit <- function(i){
        co1 <- gRbase::combnPrim(v1,i)
        co1 <- colSums(co1)
        co1 <- co1 + v[1]
        return(co1)
    }
### logits
### first element (intercept only)
    res[1,1] <- v[1]
### middle elements
    res[-c(1,nrow(res)),1] <- unlist(sapply(1:length(ch2),FUN=genLogit))
### last element (intercept + all coefficients)
    res[nrow(res),1] <- sum(v)
###
### convert logits to probabilities
    logitToProb <- function(i) exp(i)/ ( 1+exp(i) )
    res[ ,4] <- sapply(1:nrow(res), logitToProb)
### convert probabilities to Odds Ratios
    probToOR <- function(p) p/(1-p)
    res[ ,3] <- sapply(res[ ,4], probToOR)
### convert probabilities to log of Odds Ratios
    probToLnOR <- function (p) log( p/(1-p))
    res[ ,2] <- sapply(res[ ,4], probToLnOR)
###
    return(res)
}
