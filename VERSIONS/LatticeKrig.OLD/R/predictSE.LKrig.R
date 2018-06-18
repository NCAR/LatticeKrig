# fields, Tools for spatial data
# Copyright 2004-2009, Institute for Mathematics Applied to Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"predictSE.LKrig" <- function(object, xnew = NULL, 
    Znew = NULL, drop.Z = FALSE, verbose = FALSE, ...) {
    if (is.null(object$Mc)) {
        stop("need to include the sparse cholesky decompostion in LKrig object \r\nin calling LKrig set return.cholesky = TRUE")
    }
    # default is to predict at data x's
    if (is.null(xnew)) {
        xnew <- object$x
    }
    # set some local variables
    NG <- nrow(xnew)
    lambda <- object$lambda
    rho <- object$rho.MLE
    sigma2 <- lambda * rho
    weights <- object$weights
    LKinfo <- object$LKinfo
    # figure out if extra covariates should be included
    # and build fixed effects matrices
    # (ind.drift has the indices of just the spatial drift part --- e.g. the linear polynomial)
    if (is.null(Znew) & (object$nZ > 0) & (!drop.Z)) {
        Znew <- object$Z
    }
    T.matrix <- LKrig.fixed.component(object$x, Z = object$Z, 
        m = 2, distance.type = LKinfo$distance.type)
    #
    if (drop.Z | object$nZ == 0) {
        t0 <- t(fields.mkpoly(xnew, m = 2))
        Omega <- object$Omega[object$ind.drift, object$ind.drift]
    }
    else {
        t0 <- t(LKrig.fixed.component(xnew, Z = Znew, m = 2, 
            distance.type = LKinfo$distance.type))
        Omega <- object$Omega
    }
    #
    k0 <- LKrig.cov(object$x, xnew, LKinfo = LKinfo)
    wS <- diag.spam(sqrt(weights))
    PHI <- LKrig.basis(object$x, LKinfo)
    hold <- LKrig.coef(object$Mc, wS %*% LKrig.basis(object$x, 
        LKinfo), sqrt(weights) * T.matrix, sqrt(weights) * k0, 
        lambda, weights)
    # version of 'c'coefficents for usual Kriging model
    #    c.mKrig<-  weights*(k0 - T.matrix%*%hold$d.coef - PHI%*%hold$c.coef)/ lambda
    if (drop.Z | object$nZ == 0) {
        d.coef <- hold$d.coef[object$ind.drift, ]
    }
    else {
        d.coef <- hold$d.coef
    }
    # colSums used to find multiple quadratic forms  e.g.  diag(t(x) %*%A%*%x) == colSums( x* (A%*%(x)))
    temp1 <- rho * (colSums(t0 * (Omega %*% t0)) - colSums(k0 * 
        hold$c.mKrig) - 2 * colSums(t0 * d.coef))
    # find marginal variances -- trival in the stationary case!
    temp0 <- rho * LKrig.cov(xnew, LKinfo = LKinfo, marginal = TRUE)
    # Add marginal variance to part from estimate
    temp <- temp0 + temp1
    return(sqrt(temp))
    #    return(list(se=sqrt(temp), k0=k0, c.mKrig=c.mKrig, coef=hold, temp1=temp1))
}
