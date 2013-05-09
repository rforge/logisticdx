##' @name logiFactLev
##' @export
##' @title Plots factors and gives tables of factors by level
##' @description
##' For a given data frame, will check which predictors are \dfn{factors} i.e.
##' those with more with no. levels >2 and <\code{lev}.
##' \cr \cr
##' For each factor a plot is returned with \% where \emph{y}=1.
##' A (1-\code{alpha})\% confidence interval is given using Jerrreys prior.
##' \cr \cr
##' A contingency table of outcome by factor level is plotted and returned.
##' \cr \cr
##' Tukeys test is generated from \code{multcomp::glht}.
##'
##' @param df A dataframe, with predictor variables (binary, factors and
##' continuous) and a binary outcome
##' @param lev no. of \dfn{levels} which defines a factor.
##' (I.e. column in data frame is a factor if it has
##' 2< levels \eqn{\leq}{<=} \emph{lev})
##' @param alpha Significance level \eqn{\alpha}{alpha}
##' @param y Name of the outcome variable
##' @param tukey If \code{TRUE}, will add letters below plots giving Tukeys
##' all-pairwise comparison; identical letters indicate means are likely to be
##' similar by this test (at p<0.05)
##' @param round1 No. digits to which to round (for display)
##' @param width Width of screen(display device) in pixels
##' @param height Width of screen(display device) in pixels
##' @return Returns plots as detailed above, one for each factor, each on a new
##' screen.
##' \cr \cr
##' Also returns a list, with one level for each factor.
##' Each item on the list is a
##' matrix showing numbers of outcome \eqn{y=0} and \eqn{y=1} by level of factor
##' as well as \% where \eqn{y=1}, with 1-\code{alpha}\%CI given from
##' Jeffreys prior.
##' @note The letter display for Tukeys all-pairwise comparison relies on
##' \code{multcomp::cld} which can be slow for factors with e.g. >20 levels. Plots
##' work best for factors with small number of levels and can become
##' \sQuote{overcrowded} with e.g. >20 levels.
##' @seealso \code{\link{glht}} \code{\link{cld}}
##' @keywords hplot univar
##' @examples
##' set.seed(1)
##' df1 <- genLogiDf(b=1,c=1,f=3,n=20)$df
##' logiFactLev(df1,tukey=TRUE)
##'
logiFactLev <- function(df, lev=4,
                        alpha=0.05, y="y",
                        tukey=FALSE, round1=2, width=1366, height=768) {
### initialize variables (for R CMD check)
    level <- upper <- lower <- prop <- theme_text <- list1 <- NULL
###
    nr1 <- nrow(df)
### check which are binary or  continuous variables
    isFac <- function(i) isTRUE(length(unique(df[,i]))==2) |
        isTRUE(length(unique(df[,i])) > lev)
### logical vector (for columns)
    fac1 <- sapply(1:ncol(df),isFac)
### remove them
    df2 <- cbind( df[!fac1], df[y])
### ensure df contains at least one factor
    ncol1 <- ncol(df2)-1
    exp <- paste0("No factors with < ", lev, " levels")
    if (isTRUE(ncol1==0)) return(exp)
### names of factors
    fnam <- colnames(df2)[!colnames(df2)==y]
###
    plot1 <- function(j) {
### subset (per factor)
        s1 <- cbind(as.factor(df2[[j]]),df2[ncol(df2)])
        colnames(s1) <- c(colnames(df2[j]),colnames(df2[ncol(df2)]))
### make table
        t1 <- table( s1 )
        t1 <- cbind(t1, rowSums(t1))
        colnames(t1) <- c("y=0","y=1","total")
### levels of factor
        nrow1 <- length(dimnames(t1)[[1]])
        nam1 <- list(1:nrow1, c("level",
                                "y=0","y=1",
                                "total","prop","lower","upper"))
### matrix to hold CIs
        m1 <- data.frame(matrix(ncol=7, nrow=nrow1, NA, dimnames=nam1))
### levels
        m1[ ,1] <- factor(sort(unique(s1[ ,1])))
        m1[ ,2:3] <- table(s1)
        m1[ ,4] <- m1[ ,2] + m1[ ,3]
        m1[ ,5] <- round( 100*m1[,3]/m1[,4], round1)
### use alpha for confidence interval
        ci1 <- 1-alpha
### lower
        low1 <- function(i) 100*
            MKmisc::binomCI(m1[i,3],
                            m1[i,4],
                            conf.level=ci1,
                            method="jeffreys")[[2]][1]
        m1[ ,6] <- round( sapply(1:nrow1, low1), round1)
### upper
        up1 <- function(i) 100*
            MKmisc::binomCI(m1[i,3],
                            m1[i,4],
                            conf.level=ci1,
                            method="jeffreys")[[2]][2]
        m1[ ,7] <- round( sapply(1:nrow1,up1), round1)
###
        fmla <- as.formula(paste0("y ~",colnames(s1)[1]))
        f1 <- glm(fmla, family = binomial("logit"), data=s1)
        lab1 <- m1[ ,1]
        xlab1 <- "Factor level"
###
        if (isTRUE(tukey)){
### need to generate argument to linfct before using it in glht
            args <- list("Tukey")
            names(args) <- colnames(s1)[1]
            cmp1 <- do.call(multcomp::mcp, args)
            tuk <- multcomp::glht(f1, linfct=cmp1)#
### can be slow for factors with >20 levels
            c1 <- multcomp::cld(tuk, level=alpha)
### extract names and values
            df3 <- data.frame(cbind(names(c1$mcletters[[1]]),c1$mcletters[[1]]))
### convert factor back to numeric
            df3$X1 <- as.numeric(df3$X1)
### order by levels of factor
            df3 <- df3[order(df3$X1), ]
            genLab <- function(i) paste(df3[i,1],levels(df3[i,2]),sep="\n ")
            lab1 <- sapply(1:nrow(df3),genLab)
            xlab1 <- "Factor \n Letters indicate means likely to be similar \n (by Tukey's all-pair comparison with p<0.05)"
        }
###
        tit1 <- "% cases where y=1 by level of factor \n 95% CIs from Jeffreys interval"
###
     g1 <- ggplot2::ggplot(m1, ggplot2::aes(x=level, y=prop, colour=level, group=level)) +
 ggplot2::geom_point(ggplot2::aes(size=3), show_guide = FALSE) +
 ggplot2::geom_errorbar(ggplot2::aes(ymin=lower, ymax=upper), width=.1) +
 ggplot2::scale_x_discrete(xlab1, breaks=1:nrow(df3),labels=lab1 ) +
 ggplot2::scale_y_continuous("% cases \n where y=1") +
 ggplot2::scale_colour_discrete(name="Level") +
 ggplot2::opts(title=tit1,
 axis.title.x = ggplot2::theme_text(face="bold", size=15),
 axis.text.x = ggplot2::theme_text(face="bold", size=10),
 axis.title.y = ggplot2::theme_text(face="bold", size=15),
 axis.text.y = ggplot2::theme_text(face="bold", size=10),
 plot.title = ggplot2::theme_text(face="bold", size=15),
 legend.title = ggplot2::theme_text(face="bold", size=10),
 legend.text = ggplot2::theme_text(face="bold", size=10)
 )

        g2 <- gridExtra::tableGrob(m1, name="Outcome by factor")

        tit1 <- paste("Factor", colnames(df2)[j], sep=" ")
        gridExtra::grid.arrange(ggplot2::ggplotGrob(g1), g2,
                                main=grid::textGrob(tit1,
                                gp=grid::gpar(font=2, cex=2)),
                                ncol=2)
        list1 <- list(m1)
        names(list1)[1] <- colnames(df2)[j]
        return(list1)
}
    windows(record=TRUE, width=width,height=height)
    res <- sapply(1:(ncol(df2)-1), FUN=plot1)
    return(res)
}







