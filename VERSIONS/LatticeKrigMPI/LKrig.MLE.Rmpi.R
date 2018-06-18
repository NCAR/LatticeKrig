LKrig.MLE.Rmpi <- function(x, y, ..., LKinfo, par.grid = NULL, 
                  lambda.profile = TRUE, verbose = FALSE) {
  LKrig.args <- c(list(x = x, y = y), list(...))
# compare parameter grid to the LKinfo object and possibky adjust the parameter grid for
# parameter that have not been specified.
  par.grid <- LKrig.make.par.grid(par.grid = par.grid, LKinfo = LKinfo)
  if (verbose) {
    print(par.grid)
  }
  # output matrix
  NG <- length(par.grid$alpha)
  if (verbose) {
    print(NG)
  }
  if (length(par.grid$llambda) != NG) {
    stop("llambda values not same length as alpha")
  }
  if (length(par.grid$a.wght) != NG) {
    stop("a.wght values not same length as alpha")
  }
  out <- matrix(NA, nrow = NG, ncol = 9)
  optim.counts <- rep(NA, 2)
  dimnames(out) <- list(NULL, c("EffDf", "lnProfLike", "GCV", 
                                "sigma.MLE", "rho.MLE", "llambda.MLE", "lnLike", "counts value", 
                                "grad"))
  # first fit to get cholesky symbolic decomposition
  LKinfo.temp <- LKinfo
  LKinfo.temp$a.wght <- (par.grid$a.wght[[1]])
  LKinfo.temp$alpha <- (par.grid$alpha[[1]])
  LKinfo.temp$nu <- par.grid$nu[1]
  # save decomp
  MC.save <- do.call("LKrig", c(LKrig.args, list(LKinfo = LKinfo.temp, 
                                                 lambda = 1, NtrA = 0)))$MC
  
  doOneLnProfileLike<-function(k) {
    # Note each component of alpha and a.wght is also a list!
    LKinfo.temp$a.wght <- par.grid$a.wght[[k]]
    LKinfo.temp$alpha <- par.grid$alpha[[k]]
    LKinfo.temp$nu <- par.grid$nu[k]
    
    # function used in optimizer, must be in here to eval the vars in the correct scope
    temp.fn <- function(x) {
      lnLike <- do.call("LKrig", c(LKrig.args, 
                                   list(LKinfo = LKinfo.temp, 
                                        lambda = exp(x), use.cholesky = MC.save, NtrA = 0)))$lnProfileLike.FULL
    }
    
    if (lambda.profile) {
      look <- optim(par.grid$llambda[k], temp.fn, method = "BFGS", 
                    control = list(fnscale = -1, parscale = 0.1, 
                                   ndeps = 0.01, reltol = 1e-06))
      llambda.opt <- look$par
      optim.counts <- look$counts 
    }
    else {
      llambda.opt <- par.grid$llambda[k]
    }
    obj <- do.call("LKrig", c(LKrig.args, list(LKinfo = LKinfo.temp, 
                                               lambda = exp(llambda.opt), use.cholesky = MC.save, 
                                               NtrA = 20)))
    obj$llambda.opt=llambda.opt
    obj$optim.counts=optim.counts
    return(obj)
  }

  RmpiTest=FALSE;
  if (is.loaded("mpi_initialize") & (mpi.comm.size(1)>0) ) RmpiTest=TRUE;
  print( RmpiTest)
  if (RmpiTest) {
    
    cat("Parallel evalution ", fill=TRUE) 
    par.grid$llambda[is.na(par.grid$llambda)]=par.grid$llambda[1] # can't do sequential starting points
    
    # Here's the operation locally:
    #    objAll<-lapply(1:NG,doOneLnProfileLike)
    
    # load up the slave environments
    mpi.bcast.cmd(library(LatticeKrig))
    mpi.bcast.Robj2slave(par.grid)
    mpi.bcast.Robj2slave(lambda.profile)
    mpi.bcast.Robj2slave(LKinfo.temp)
    mpi.bcast.Robj2slave(MC.save)
    
    # do all the calculations
    objAll=mpi.applyLB(1:NG,doOneLnProfileLike)
  }
  
  lnProfileLike.max <- -Inf
  for (k in 1:NG) {
    # if starting value is missing use the optimum from previous fit
    # this only really makes sense if other parameters have some sort of continuity from k-1 to k.
    
    if (RmpiTest) {
      obj<-objAll[[k]]
    }
    else
    {
      if (is.na(par.grid$llambda[k]) & (k != 1)) {
        par.grid$llambda[k] <- llambda.opt
      }      
      obj<-doOneLnProfileLike(k);
    }
    llambda.opt=obj$llambda.opt
    # compare to current largest likelihood and update the LKinfo list if bigger.
    if (obj$lnProfileLike.FULL > lnProfileLike.max) {
      lnProfileLike.max <- obj$lnProfileLike.FULL
      LKinfo.MLE <- obj$LKinfo
    }
    #
    out[k, 1] <- obj$trA.est
    out[k, 2] <- obj$lnProfileLike.FULL
    out[k, 3] <- obj$GCV
    out[k, 4] <- obj$sigma.MLE.FULL
    out[k, 5] <- obj$rho.MLE.FULL
    out[k, 6] <- llambda.opt
    out[k, 7] <- obj$lnLike.FULL
    out[k, 8:9] <- obj$optim.counts
    if (verbose) {
      cat("  ", k, "eff.df:", out[k, 1], "lnProfile like:", 
          out[k, 2], "llambda:", out[k, 6], "lnlike:", 
          out[k, 7], fill = TRUE)
    }
  } # for (k in 1:NG)
  
  return(list(summary = out, par.grid = par.grid, LKinfo = LKinfo, 
              LKinfo.MLE = LKinfo.MLE, call = match.call()))
}
