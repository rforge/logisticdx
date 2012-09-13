genLogiDf <-
function(b=2, f=0, c=1, n=10, rb=0.5, rf=0.3, rc=0.8, ry=0.5, asFactor=FALSE, timeout=5) {
###require(speedglm)
### for debugging:
###    b=0; f=2; c=0; n=5000; rb=0.5; rf=0.1; rc=0.8; ry=0.5; asFactor=TRUE; timeout=10

    if ( (rb|rf|rc|ry)>1 ) return("Ratios should be <1")
    bcf <- c(b,c,f)
    if (sum(bcf) <= 0) return("Need at least one predictor")
    if (sum(bcf)>=n) return("Need n to be larger for this no. predictors")

setTimeLimit(elapsed=timeout, transient=TRUE) # prevent taking more than 10secs
    repeat{
    y1 <- sample(x=c(0,1),size=n,replace=TRUE,prob=c(ry,1-ry))
    b1 <- sample(x=c(0,1),size=b*n,replace=TRUE,prob=c(rb,1-rb))
    b1 <- tryCatch( matrix(b1, ncol=(b)), error=function(e) NA)
    tryCatch( colnames(b1) <- rep("b",ncol(b1)), error=function(e) NA)
    c1 <- sample(x=seq(1,(rc*n)),size=c*n,replace=TRUE)
    c1 <- tryCatch( matrix(c1, ncol=c), error=function(e) NA)
    tryCatch( colnames(c1) <- rep("c",ncol(c1)), error=function(e) NA)
    f1 <- sample(x=seq(1,(rf*n)),size=f*n,replace=TRUE )
    if (length(unique(f1))==2) { # rename as 0,1 if only 2x elements
        tryCatch( f1[f1==min(f1)] <- 0,  error=function(e) NA)
        tryCatch( f1[f1==max(f1)] <- 1,  error=function(e) NA)
        }
    f1 <- tryCatch( matrix(f1, ncol=f), error=function(e) NA)
    ifelse( length(unique(as.vector(f1)))==2, # rename as b if only 2 elements
           tryCatch( colnames(f1) <- rep("b",ncol(f1)), error=function(e) NA),
           tryCatch( colnames(f1) <- rep("f",ncol(f1)), error=function(e) NA)
           )
    w1 <- cbind(b1,f1,c1,y1)
    w1 <- w1[,colMeans(is.na(w1)) == 0] # remove columns with NA
    storage.mode(w1) <- "integer" # duplicated function doesn't work with double
    if (anyDuplicated(w1, MARGIN = 2) == FALSE) {break} # check no duplicate columns
    }
    w1 <- data.frame(w1)
    if (asFactor==TRUE) w1[,c(grep("f",colnames(w1)))] <- lapply( w1[,c(grep("f",colnames(w1)))], as.factor)
    xnam <- paste("x", 1:(ncol(w1)-1), sep="")
    colnames(w1) <- c(xnam,"y")

    fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")  ))
    f1 <- glm (formula=fmla, family=binomial(), data=w1) # All predictors
###f1 <- speedglm::speedglm (formula=fmla, family=binomial(), data=w1) # All predictors

res <- list(df=w1, model=f1)
return(res)
}
