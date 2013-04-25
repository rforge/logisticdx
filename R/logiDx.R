logiDx <-
function(model, width=TRUE){

    if (class(model)[1] !="glm") return("logiDx only applies to objects of class glm")


v <- model$coefficients
l1 <- length(v)-1 # no. coefficients, excluding intercept

w1 <- model$data

### make design/ covariate pattern matrix; all combinations of predictors

fun1 <- function(j) paste("length(unique(w1[,",j,"]))",sep="") # creates expression for each unique subsets
e1t <- sapply(1:l1,fun1)
e2 <- paste(e1t, collapse="*")
nr1 <- eval(parse(text=e2))  # no rows = no covariate patterns
nc1 <- l1 + 11 # no ools; 11x additional measurements
m1 <- matrix(NA,ncol=nc1, nrow=nr1) # matrix to fill

fun1 <- function(j) paste("unique(w1[,",j,"])",sep="") # creates expression for each unique subsets
e1t <- sapply(1:l1,fun1)
e1 <- paste(e1t, collapse=" , ")
e1 <- paste("expand.grid(",e1,")", sep="")
e1 <- parse(text=e1) # convert character string to expression
m1[(1:nr1),(1:l1)] <- as.matrix(eval(e1))
colnames(m1) <- 1:dim(m1)[2] # arbitrary staring column names
colnames(m1)[1:(l1)] <- names(v)[2:length(v)]

### no observations/co-variate pattern
fun1 <- function(j) paste("w1[names(v)[2:length(v)][",j,"]]==m1[i,",j,"]", sep="") # creates expression for each unique subset
e1t <- sapply(1:l1, fun1)
e1 <- paste(e1t, collapse=" & ")
e1 <- parse(text=e1) # convert character string to expression
fun1 <- function(i) nrow( subset(w1, eval(e1)))
m1[,l1+1] <- sapply(1:nr1, fun1) # obs is no observations fitting criteria

### probability (predicted) each co-variate pattern
fun1 <- function(j) paste("v[",j+1,"] * m1[i,",j,"]", sep="") # creates expression for each unique subset
e1t <- sapply(1:l1, fun1)
e1 <- paste(e1t, collapse=" + ")
e1 <- paste("v[1]",e1, sep=" + ")
e1 <- parse(text=e1) # convert character string to expression
fun1 <- function(i) round( exp(eval(e1))/ (1+exp(eval(e1))), 3) # eval(e1) gives logit for this covariate pattern
m1[,l1+2] <- sapply(1:nr1,fun1) # pi is predicted probability for each pattern

### yhat =  predicted no. outcome=1 in each pattern
m1[,l1+3] <- m1[,l1+1]*m1[,l1+2] # yhat = obs * pi

###acutal no.
e1 <- vector(length = l1, mode="character")
###fun1 <- function(j) paste("w1[xnam[",j,"]]==m1[i,",j,"]", sep="") # creates expression for each unique subsets
fun1 <- function(j) paste("w1[names(v)[2:length(v)][",j,"]]==m1[i,",j,"]", sep="") # creates expression for each unique subsets
e1 <- sapply(1:l1, fun1)
e1 <- paste(e1, collapse=" & ")
e1 <- parse(text=e1);e1 # convert character string to expression
### actual no observations fitting these criteria:
fun1 <- function(i) {
	ss <- subset(w1,eval(e1),select="y")
	return( sum(ss$y==1) )
	}
m1[,l1+4] <- sapply(1:nr1,FUN=fun1)

colnames(m1)[(l1+1):(l1+4)] <- c("obs","pi","yhat","y")

###--------------------------------
### make hat matrix
d1 <- cbind(rep(1,nrow(m1)), m1[,(1:l1)] ) # add intercept term to make Design matrix

fun1 <- function(i) m1[i,"obs"] * m1[i,"pi"] * (1-(m1[i,"pi"]))
v <- sapply(1:nr1,FUN=fun1) #no observations fitting criteria
v1 <- diag(v) # make diagonal matrix
v1s <- sqrt(v1)
H1 <- v1s %*% d1 %*% solve(t(d1) %*% v1 %*% d1) %*% t(d1) %*% v1s # short form
m1[,l1+5] <- round( diag(H1), 3) # hat diagonals = leverage
rm (v, v1, v1s)

### deviance residual for covariate pattern
fun1 <- function(j){
	if (m1[j,"y"] ==0){ # y=0 for this covariate pattern
		d1 <- log( (1-m1[j,"pi"]) )
		dev <- -sqrt( 2 * m1[j,"obs"] * abs(d1) )
		res <- round(dev, 3)
		return(res)
		} else if (m1[j,"y"] == m1[j,"obs"]){
			d1 <- log( (m1[j,"pi"]) )
			dev <- sqrt( 2 * m1[j,"obs"] * abs(d1) )
			res <- round(dev, 3)
			return(res)
			} else {
	d1 <- m1[j,"y"] / (m1[j,"yhat"])
	d2 <- m1[j,"y"] * log(d1)
	d3 <- ( m1[j,"obs"] - m1[j,"y"] ) / ( m1[j,"obs"] * (1 -m1[j,"pi"]) );d3
	d4 <- (m1[j,"obs"]-m1[j,"y"] )* log(d3)
	d5 <- sqrt( 2*(d2 + d4) )
	s1 <- sign(m1[j,"y"]- m1[j,"yhat"]) # 1 if +ve, 0 if -ve
	dev <- s1*d5; dev
	res <- round(dev, 3)
	return(res)
	}
}
m1[,l1+6] <- sapply(1:nr1,fun1) # devR

### Pearson residual
fun1 <- function(j){
	if (m1[j,"y"] ==0){ # y=0 for this covariate pattern
		Pr1 <- sqrt( m1[j,"pi"] / (1-m1[j,"pi"]))
		Pr2 <- -sqrt (m1[j,"obs"])
		res <- round( Pr1 * Pr2, 3);res
		return(res)
		} else {
	Pr1 <- m1[j,"y"] - m1[j,"yhat"]
	Pr2 <- sqrt(   m1[j,"yhat"] * ( 1-(m1[j,"pi"]) )   )
	res <- round( Pr1/Pr2, 3);res
	return(res)
	}
}
m1[,l1+7] <- sapply(1:nr1, FUN=fun1) # PeR

colnames(m1)[(l1+5):(l1+7)] <- c("lev","devR","PeR")

### Standardized Pearson residuals
fun1 <- function(i) round ( m1[i,"PeR"] / (sqrt (1-m1[i,"lev"]) ), 3)
m1[,l1+8] <- sapply(1:nr1,FUN=fun1) #sPeR

colnames(m1)[(l1+8)] <- "sPeR"

### standardized difference in B (maximum likelihood coefficients for paramaters) without this pattern (standardized by covariance matrix of B); should be <1 if little influence on model
fun1 <- function(i) round ( (m1[i,"sPeR"])^2 * m1[i,"lev"] / (1-m1[i,"lev"]), 3)
m1[,l1+9] <- sapply(1:nr1, FUN=fun1) #dBhat

### decrease in Pearson chi-square without this pattern, should be <4 if little influence on model
fun1 <- function(i) round ( (m1[i,"sPeR"])^2, 3)
m1[,l1+10] <- sapply(1:nr1, FUN=fun1) #dXsq

### decrease in Pearson chi-square without this pattern, should be <4 if little influence on model
fun1 <- function(i) round ( (m1[i,"devR"])^2 / (1-m1[i,"lev"]) , 3)
m1[,l1+11] <- sapply(1:nr1, FUN=fun1) #dDev

colnames(m1)[(l1+9):(l1+11)] <- c("dBhat","dXsq","dDev")

if(isTRUE(width)) options(width= ((l1+11)*6) ) # allows all to be seen together

res <- list(dxMatrix=m1)
###class(res) <- c("logisticDx","matrix")
return(res)
}
