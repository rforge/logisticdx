logiFactLev <-
function(df, r=0.4, y="y", tukey=FALSE, width=1366, height=768) {
### require(multcomp)
### require(gridExtra)
### require(ggplot2)
### require(MKmisc)
### for debugging:
###    r=0.45; y="y"; tukey=TRUE; width=1366; height=768

    nr1 <- nrow(df)

### check which are binary or  continuous variables

    fun1 <- function(i) isTRUE(length(unique(df[,i]))==2) | isTRUE(length(unique(df[,i]))/nr1 > r)
    w2 <- sapply(1:ncol(df),fun1);w2
    w2 <- cbind( df[!w2], df[y]) # remove them

    ncol1 <- ncol(w2)-1
    exp <- paste("No factors with < ", round(r*nr1,0), " levels", sep="")
    if (isTRUE(ncol1==0)) return(exp)

    fnam <- paste("f", 1:ncol1, sep="")
    colnames(w2) <- c(fnam, "y")

fun1 <- function(j) {

    s1 <- w2[,c(j,ncol(w2))]            # subset
    colnames(s1) <- c("f1","y")

    t1 <- table( s1[,c(1,2)] )
    t1 <- cbind(t1, rowSums(t1))
    colnames(t1) <- c("y=0","y=1","total")

    nrow1 <- length(dimnames(t1)[[1]])  # levels of factor
    nam1 <- list(1:nrow1, c("level","y=0","y=1","total","prop","lower","upper"))
    m1 <- matrix(ncol=7,nrow=nrow1, NA, dimnames=nam1) # matrix to hold CIs
    m1[,1] <- as.numeric(dimnames(t1)[[1]]) # levels
    m1[,2:4] <- t1                          # table above
    m1[,5] <- round( 100*m1[,3]/m1[,4], 2)

    fun1 <- function(i) 100* MKmisc::binomCI(m1[i,3], m1[i,4], conf.level=0.95, method="jeffreys")[[2]][1]
    m1[,6] <- round( sapply(1:nrow1,fun1), 2)
    fun1 <- function(i) 100* MKmisc::binomCI(m1[i,3], m1[i,4], conf.level=0.95, method="jeffreys")[[2]][2]
    m1[,7] <- round( sapply(1:nrow1,fun1), 2)

    list1 <- list(m1)
    names(list1)[1] <- colnames(w2)[j]

    m1 <- data.frame(m1)                # change to data frame to plot
    m1$level <- factor(m1$level)

    s1[,1] <- factor(s1[,1])
    fmla <- as.formula(y~f1)
    f1 <- glm (fmla, family = binomial("logit"), data=s1)

    lab1 <- m1[,1]
    xlab1 <- "Factor level"
    if (isTRUE(tukey)){
        tuk <- multcomp::glht(f1,linfct=multcomp::mcp(f1="Tukey"))
        c1 <- multcomp::cld(tuk, level=.05) # can be slow for factors with >20 levels
        df1 <- data.frame(cbind(names(c1$mcletters[[1]]),c1$mcletters[[1]]))
        df1$X1 <- as.numeric(df1$X1)
        df1 <- df1[order(df1$X1),]      # order by levels of factor
        fun1 <- function(i) paste(df1[i,1],levels(df1[i,2]),sep="\n ")
        lab1 <- sapply(1:nrow(df1),fun1)
        xlab1 <- "Factor \n Letters indicate means likely to be similar \n (by Tukey's all-pair comparison with p<0.05)"
    }

    tit1 <- "% cases where y=1 by level of factor \n 95% CIs from Jeffreys interval"

    g1 <- ggplot2::ggplot(m1, ggplot2::aes(x=level, y=prop, colour=level, group=level)) +
 ggplot2::geom_point(ggplot2::aes(size=3), show_guide = FALSE) +
 ggplot2::geom_errorbar(ggplot2::aes(ymin=lower, ymax=upper), width=.1) +
 ggplot2::scale_x_discrete(xlab1, breaks=1:nrow(df1),labels=lab1 ) +
 ggplot2::scale_y_continuous("% cases \n where y=1") +
 ggplot2::scale_colour_discrete(name="Level") +
 ggplot2::opts(title=tit1,
 axis.title.x = theme_text(face="bold", size=15),
 axis.text.x = theme_text(face="bold", size=10),
 axis.title.y = theme_text(face="bold", size=15),
 axis.text.y = theme_text(face="bold", size=10),
 plot.title = theme_text(face="bold", size=15),
 legend.title = theme_text(face="bold", size=10),
 legend.text = theme_text(face="bold", size=10)
 )

    g2 <- gridExtra::tableGrob(m1, name="Outcome by factor")

    tit1 <- paste("Factor", colnames(w2)[j], sep=" ")
    gridExtra::grid.arrange(g1, g2, main=grid::textGrob(tit1, gp=grid::gpar(font=2, cex=2)), ncol = 2)
return(list1)
}

    windows(record=TRUE, width=width,height=height)
    res <- sapply(1:(ncol(w2)-1),FUN=fun1)
    return(res)
}
