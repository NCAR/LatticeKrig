LKrigFindLambda <- function(x, y, ..., LKinfo,
                            use.cholesky=NULL, 
                            verbose = FALSE,
                            lambda.profile=TRUE,
                            lowerBoundLogLambda=-16) {
    LKrigArgs <- c(list(x = x, y = y),
                   list(...),
                   list( LKinfo=LKinfo, NtrA=ifelse( lambda.profile, 0, 20) ))
    out <- rep(NA, 9)
    names( out) <-  c("EffDf", "lnProfLike", "GCV", 
        "sigma.MLE", "rho.MLE", "llambda.opt", "lnLike", "counts value", "grad")
    capture.evaluations <- matrix(NA, ncol = 4, nrow = 1, 
                dimnames = list(NULL, c("lambda", "rho.MLE", "sigma.MLE", "lnProfileLike.FULL")))
    optim.counts<-  NA
    llambda.start<- log( LKinfo$lambda )
# first fit to get cholesky symbolic decomposition  and wPHI matrix
# Note that if use.cholesky != NULL then the symbolic decomposition information is
# used from the passed object.
    if( is.na(llambda.start)){ llambda.start<- -1 }
    LKrigObject <- do.call("LKrig", c(
                             LKrigArgs,
                             list( lambda = exp( llambda.start),
                                   use.cholesky=use.cholesky,
                                   return.cholesky = TRUE,
                                   return.wPHI=TRUE,
                                   verbose=FALSE)))
    capture.evaluations[1,] <-  unlist(LKrigObject [c("lambda", "rho.MLE", "sigma.MLE", "lnProfileLike.FULL")])
    llambda.opt<- llambda.start
    Mc.save<- LKrigObject$Mc
    wPHI.save<- LKrigObject$wPHI
#    
##### in most cases now optimze likelihood over log lambda    
   if( lambda.profile){
#  temporary function used in optimizer
#   
# objective function    
    temp.fn <- function(x) {
                             hold <- do.call("LKrig", c(LKrigArgs,
                                      lambda = exp(x), use.cholesky = Mc.save, wPHI= wPHI.save, verbose=FALSE)
                                             )[c("lambda", "rho.MLE", "sigma.MLE", "lnProfileLike.FULL")]
            temp.eval <- get("capture.evaluations")
            assign("capture.evaluations", rbind(temp.eval, unlist(hold)), envir = capture.env)
                             return(hold$lnProfileLike.FULL)}
#    
# the try wrapper captures case when optim fails.   
            capture.env <- environment()
            look<- try(optimize( temp.fn, interval = llambda.start+c(-8,5), maximum= TRUE, tol=1e-3))
            evalSummary <- !(class( look)== "try-error")
            llambda.opt <- look$maximum
            optim.counts <- nrow( capture.evaluations) + 2
            LKrigArgs$NtrA<- 20
            LKinfo$lambda<- exp(llambda.opt)
            LKrigObject <- do.call("LKrig", c(LKrigArgs, use.cholesky = Mc.save,
                                             wPHI= wPHI.save, lambda = exp(llambda.opt)))           
  }
###### end optimze block    
# save summary results from this set of parameters.
   out[ 1] <- LKrigObject$trA.est
   out[ 2] <- LKrigObject$lnProfileLike.FULL
   out[ 3] <- LKrigObject$GCV
   out[ 4] <- LKrigObject$sigma.MLE.FULL
   out[ 5] <- LKrigObject$rho.MLE.FULL
   out[ 6] <- llambda.opt
   out[ 7] <- LKrigObject$lnLike.FULL
   out[ 8] <- optim.counts
   out[ 9] <- NA
   return(list(summary = out,
                LKinfo = LKinfo,
                llambda.start= llambda.start,
                lambda.MLE =  exp(llambda.opt),
                lnLike.eval=capture.evaluations,
                call = match.call(),
                Mc=Mc.save)
          )
}

