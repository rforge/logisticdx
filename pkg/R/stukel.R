stukel <-
function(object, alternative = c("both", "alpha1", "alpha2")) {
  DNAME <- deparse(substitute(object))
  METHOD <- "Stukel's test of the logistic link"
  alternative <- match.arg(alternative)
  eta <- predict(object, type = "link")
  etasq <- 0.5 * eta * eta
  etapos <- eta > 0
  dv <- matrix(0, nrow = length(eta), ncol = 2)
  dv[etapos,1] <- etasq[etapos]
  dv[!etapos,2] <- - etasq[!etapos]
  colnames(dv) <- c("z1","z2")
  oinfo <- vcov(object)
  oX <- qr.X(object$qr)
  ImH <- - oX %*% oinfo %*% t(oX)
  diag(ImH) <- 1 + diag(ImH)
  wdv <- sqrt(object$weights) * dv
  qmat <- t(wdv) %*% ImH %*% wdv
  sc <- apply(dv * (object$weights * residuals(object, "working")), 2, sum)
  allstat <- c(sc * sc / diag(qmat), sc %*% solve(qmat) %*% sc)
  names(allstat) <- c("alpha1", "alpha2", "both")
  allpar <- c(1,1,2)
  names(allpar) <- names(allstat)
  allpval <- pchisq(allstat, allpar, lower.tail=FALSE)
  STATISTIC <- allstat[alternative]
  PARAMETER <- allpar[alternative]
  names(PARAMETER) <- "df"
  PVAL <- allpval[alternative]
  names(allpar) <- rep("df", 3)
  structure(list(statistic = STATISTIC,
                 parameter = PARAMETER,
                 p.value = PVAL,
                 alternative = alternative,
                 method = METHOD, data.name = DNAME,
                 allstat = allstat, allpar = allpar, allpval = allpval
                 ),
            class = "htest")
}
