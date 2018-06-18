LKrig.sim.conditional <- function(obj, M = 1, x.grid = NULL, 
    grid.list = NULL, nx = 80, ny = 80, ..., Z.grid = NULL) {
    # generate grid if not specified
    if (is.null(x.grid)) {
        if (is.null(grid.list)) {
            grid.list <- fields.x.to.grid(obj$x, nx = nx, ny = ny)
        }
        x.grid <- make.surface.grid(grid.list)
    }
    ghat <- predict(obj, xnew=x.grid, Znew = Z.grid)
    # NOTE the name x.grid may be misleading because it just needs to a column matrix of
    # locations. It need not follow any regualr pattern.
    # now generate the error surfaces
    # begin block
    g.conditional.draw <-    matrix(NA, ncol = M, nrow = nrow(x.grid))
    d.coef.draw<- matrix(NA, ncol = M, nrow = length( obj$d.coef) )
    set.seed(122)
    N <- nrow(obj$x)
    # complete set of locations to evaluate the field must include the observations too
    X.full <- rbind(obj$x, x.grid)
    for (k in 1:M) {
        cat(k, " ")
        g.full <- LKrig.sim(X.full, LKinfo = obj$LKinfo) * sqrt(obj$rho.MLE)
        g.unconditional.data <- g.full[1:N]
        g.unconditional.grid <- g.full[-(1:N)]
        # generate a synthetic data set with fixed part set to zero.
        y.synthetic.data <- g.unconditional.data + obj$sigma.MLE * 
            rnorm(N)/obj$weights
        # use LKrig to find the predictions for the xgrid locations
        # NOTE that LKrig will still estimate the fixed part.
        obj.fit.synthetic <- LKrig(obj$x, y.synthetic.data, LKinfo = obj$LKinfo, 
            lambda = obj$lambda, Z = obj$Z, weights = obj$weights, 
            ...)
        d.coef.draw[,k] <- obj.fit.synthetic$d.coef
        #
        error <- g.unconditional.grid - predict(obj.fit.synthetic, 
                            x.grid, drop.Z = is.null(Z.grid), Znew = Z.grid)
        # this result has the right covariance but zero mean add the conditional mean to it.
        g.conditional.draw[, k] <- ghat + error
    }
    #
    return(list(x.grid = x.grid, ghat = ghat, g.draw = g.conditional.draw,
                           d.coef.draw= d.coef.draw))
}
