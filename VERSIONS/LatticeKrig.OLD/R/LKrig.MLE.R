LKrig.MLE <- function(x, y, ..., LKinfo, par.grid = NULL, 
    lambda.profile = TRUE, verbose = FALSE, lowerBoundLogLambda=-16) {
    LKrig.args <- c(list(x = x, y = y), list(...))
    # at this point LKinfo has the correct value for the number of multiresolution levels
    par.grid <- LKrig.make.par.grid(par.grid = par.grid, LKinfo = LKinfo)
    # output matrix
    NG <- length(par.grid$alpha)
    out <- matrix(NA, nrow = NG, ncol = 9,
                  dimnames = list(NULL, c("EffDf", "lnProfLike", "GCV", 
                  "sigma.MLE", "rho.MLE", "llambda.MLE", "lnLike", "counts value", "grad")))
    optim.counts <- rep(NA, 2)
    # evaluate parameters but do an optimzation over lambda
    lnProfileLike.max <- -1e+20
    lnLike.eval<- NULL
    llambda.opt<- NA
    for (k in 1:NG) {
    # if starting value is missing use the optimum from previous fit
    # this only really makes sense if other parameters have some sort of continuity from k-1 to k.
        llambda.start<- ifelse (is.na( par.grid$llambda[k]), llambda.opt,  par.grid$llambda[k] )  
     # first fit to get cholesky symbolic decomposition
        LKinfo.temp<- LKinfoUpdate( LKinfo, a.wght = (par.grid$a.wght[[k]]),
                                    alpha =  (par.grid$alpha[[k]]),
                                    nu = par.grid$nu[k], lambda=exp(llambda.start) )
    # for first pass save symbolic decomposition for M
        if( k ==1 ){ use.cholesky <- NULL}
    #    
        obj <- do.call("LKrigFindLambda",c(
                                       LKrig.args,
                                       list(LKinfo = LKinfo.temp,lambda.profile=lambda.profile)
                                           )
                                      )
     # Note: if lambda.profile == FALSE the starting lambda is passed through as the "optimal" one
        llambda.opt<- obj$summary["llambda.opt"]
     # compare to current largest likelihood and update the LKinfo.MLE list if bigger.
        if ( obj$summary["lnProfLike"] > lnProfileLike.max ) {
              lnProfileLike.max <- obj$summary["lnProfLike"] 
              LKinfo.MLE <- obj$LKinfo
              lambda.MLE<- exp( llambda.opt)
        }
     # save summary results from this set of parameters.
       out[k, ] <- obj$summary
       lnLike.eval <- c( lnLike.eval, list(obj$lnLike.eval))
     #   
       if (verbose) { print( c(k,out[k,])) }
    }
    return(list(
                summary = out,
                par.grid = par.grid,
                LKinfo = LKinfo, 
                LKinfo.MLE = LKinfo.MLE,
                lnLike.eval = lnLike.eval,
                lambda.MLE = lambda.MLE,
                call = match.call() )
           )
}

