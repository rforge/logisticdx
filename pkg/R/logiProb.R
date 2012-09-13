logiProb <-
function(model){

### (optional) require(gRbase)
### for gRbase::combnPrim instead of combn (works faster)
### the following dependencies may be necessary, install as follows:
### source("http://bioconductor.org/biocLite.R")
### biocLite("graph")
### biocLite("RBGL")

    if (class(model)[1] !="glm") return("LogiProb only applies to objects of class glm")

### useful for debugging:
### nCk <- function(n,k) gamma(n+1)/ ( gamma(k+1)*gamma((n-k)+1) )
### len1 <- nCk(le1,i) # no combinations to be returned

    v <- model$coefficients
    names(v)[1] <- "Int~"               # short for intercept
    v1 <- v[-1]                         # remove intercept
    le1 <- length(v1) # used to generate combinations without intercept
    ch1 <- choose(le1,0:le1)  # binomial expansion; give no. of logits
    len1 <- sum(ch1)          # no. patterns by combination
                                        # remove first and last values from ch1:
    ch2 <- ch1[-1]
    length(ch2) <- length(ch2)-1

    na1 <- rep("",len1)                 # holds names
    na1[1] <- names(v)[1]
    fun1 <- function(i){
        co <- combn(names(v1),i)
        fun2 <- function(j) paste(co[,j],collapse="+")
        col1 <- 1:ncol(co)
        co <- sapply(col1, FUN=fun2, USE.NAMES=FALSE)
        fun3 <- function(k) paste (names(v)[1], co[k], sep="+")
        co <- sapply(col1, FUN=fun3)
        return(co)
    }
    na1[2:(sum(ch1)-1)] <- unlist(sapply(1:length(ch2),FUN=fun1))
    na1[len1] <- paste(names(v), collapse="+")

    lo <- rep(NA,len1)                  # holds logits
    lo[1] <- v[1]
    fun1 <- function(i){
        co <- combn(v1,i)
        co <- colSums(co)
        co <- co + v[1]
        return(co)
    }
    lo[2:(sum(ch1)-1)] <- unlist(sapply(1:length(ch2),FUN=fun1))
    lo[len1] <- sum(v)

    fun1 <- function(i) exp(i)/ ( 1+exp(i) ) # convert to probabilities
    pr <- sapply(lo, fun1)

    fun1 <- function(p) p/(1-p) # convert probabilities to Odds Ratios
    or <- sapply(pr, fun1)

    fun1 <- function (p) round (log( p/(1-p)), 5)
    lor <- sapply(pr, fun1)

    res <- cbind(lo,lor,or,pr)
    rownames(res) <- na1
    colnames(res) <- c("logit","ln(OR)","OR","p")
    return (as.data.frame(res))
}
