##' @name genLogiDf
##' @export
##' @title Generate logistic data frame and model
##' @description
##' Generates a data frame with a binary outcome, and a logistic model to
##' describe it. Model is fitted with \code{glm()}.
##'
##' @param b \dfn{binomial predictors}, the number of predictors which are
##' binary, i.e. limited to 0 or 1
##' @param f \dfn{factors}, the number of predictors which are factors
##' @param c \dfn{continuous predictors}, the number of predictors which are
##' continuous
##' @param n number of observations in the data frame
##' @param nf the no. of levels in a factor
##' @param rb \dfn{ratio for binomnial predictors} the ratio of 1s to total
##' observations for the binomial predictors e.g. if \code{rb=0.3},
##'  30\% will be 1s, 70\% will be 0s
##' @param rc \dfn{ratio for continuous variables} the ratio of levels of
##' continuous variables to the total number of observations \dfn{n} e.g. if
##' \code{rc=0.8} and \code{n=100}, it will be in the range 1-80
##' @param ry \dfn{ratio for y} the ratio of 1s to total observations for the
##' binomial predictors e.g. if \code{ry=0.5},
##' 50\% will be 1s, 50\% will be 0s
##' @param asFactor if \code{TRUE}, predictors given as factors
##' will be converted to factors in data frame before model is fit
##' @param timelim function will timeout after \code{timelim} secs
##' @param speedglm return fitted model with \code{speedglm} instead of
##' \code{glm}
##' @return A list with the following values:
##'  \item{df}{data frame with predictors
##' (labelled \eqn{x1,x2, ..., xn}) and outcome
##'   (y), with \emph{n} rows (observations)}
##'  \item{model}{model fit with
##'   \code{ glm(family=binomial())} or \code{speedglm(family=binomial())} }
##' @note Using \code{asFactor=TRUE} with factors which have a large number of
##' levels (e.g. \code{nf >30}) on large datasets (e.g. \eqn{n >1000}) can cause
##' fitting to be excessively slow.
##'
##' @keywords datagen
##' @examples
##' set.seed(1)
##' genLogiDf()
##' genLogiDf(b=0,c=2,n=100,rc=0.7)
##' genLogiDf(b=1,c=0,n=1000)
##'
genLogiDf <- function(f=0, b=2, c=1, n=20,
                      rb=0.5, nf=3, rc=0.8, ry=0.5,
                      asFactor=FALSE, timelim=5, speedglm=FALSE) {
    if ( (rb|rc|ry)>1 ) stop("Ratios should be <1")
    if (nf<=2) stop("Factors should have at least 3 levels")
    bcf <- c(b,c,f)
    if (sum(bcf) <= 0) stop("Need at least one predictor")
    if (sum(bcf)>=n) stop("Need n to be larger for this no. predictors")
### prevent taking more than (timelimit)
    setTimeLimit(elapsed=timelim, transient=TRUE)
### define frame to hold values
    df1 <- as.data.frame(matrix(0L, ncol=(sum(bcf)+1),nrow=n))
    xnam <- paste("x", 1:(ncol(df1)-1), sep="")
    colnames(df1) <- c(xnam,"y")
### repeat until no dupicated columns
    repeat{
        if(!f==0) df1[,1:f] <- sample(x=seq(1,nf),size=f*n,replace=TRUE )
        if(!b==0) df1[,(f+1):(f+b)] <- sample(x=c(0L,1L), size=b*n, replace=TRUE, prob=c(rb,1-rb) )
        if(!c==0) df1[,(b+f+1):(b+f+c)] <- sample(x=seq(1,(rc*n)),size=c*n,replace=TRUE)
        df1[,ncol(df1)] <- sample(x=c(0L,1L), size=n, replace=TRUE, prob=c(ry, 1-ry) )
### check no duplicate columns (prevent overfitting)
        if (!all(duplicated(df1))){break}
    }
### convert to factor
    if (!f==0 & asFactor==TRUE) df1[,1:f] <- do.call(as.factor, list(df1[,1:f] ))
### generate formula
    fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")  ))
### use speeedglm?
    if(speedglm){
        f1 <- speedglm:::speedglm(formula=fmla, family=binomial(), data=df1)
    } else {
        f1 <- glm (formula=fmla, family=binomial("logit"), data=df1)
    }
    res <- list(df=df1,
                model=f1)
    return(res)
}
