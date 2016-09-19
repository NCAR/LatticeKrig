
# ####################################################################################
LKrigFindLambdaAwght <- function(x, y, ...,  LKinfo,
                                 use.cholesky=NULL, 
                                 lowerBoundLogLambda =-16,
                                 upperBoundLogLambda = 4,
                                 lowerBoundOmega = -5,
                                 upperBoundOmega =  2,
                                 verbose=FALSE) {
  require(stats)
  
  # For rectangle omega = log(kappa) = log(sqrt(Awght-4))
  # but will change with other models. 
  # Parts of the LKrig call that will be fixed.  (except updates to LKinfo)                             
  LKrigArgs <- c(list(x = x, y = y), list( ...),
                 list( LKinfo=LKinfo, NtrA= 0  ))
  if( verbose){
    cat( "LKrigFindLambdaAwght: Set of LKrigArgs before first call:", names(LKrigArgs ), fill=TRUE)
  }
  # Set up initial values of Awght and omega
  Awght.init <- as.numeric(LKrigArgs$LKinfo$a.wght[1])
  omega.start <- Awght2Omega( Awght.init, LKinfo)
  
  # Set up initial values of lambda and log lambda
  lambda.start <- LKrigArgs$LKinfo$lambda 
 
  if( is.na(lambda.start) ){ 
    llambda.start <- -1 
  } 
  else{
    llambda.start<- log( lambda.start)
    }
# 
  if ( (llambda.start < lowerBoundLogLambda) || 
       (llambda.start > upperBoundLogLambda) || 
       (is.na(llambda.start))
      ){
    stop("Given lambda value is out of bounds.")
  }  
#
  if(verbose){
    cat("LKrigFindLambdaAwght: llambda.start:",  llambda.start, "a.wght.start:", Awght.init, fill=TRUE)
  }
  a.wghtTemp<- Omega2Awght(omega.start, LKinfo)
  lambdaTemp<- exp( llambda.start)
  # initial call to likelihood and also to get symbolic decomposition of 
  # the "M" matrix k-- sparsity pattern does not change as lambda, awght are varied.
  LKrigObject <- do.call("LKrig", c(
    LKrigArgs,
    list( 
      use.cholesky = use.cholesky, 
      return.cholesky = TRUE,
      return.wXandwU = TRUE,
      lambda = lambdaTemp,
      a.wght = a.wghtTemp,
      verbose = FALSE)))
# Update the LKrigArgs with cholesky decomp and wU  
 LKrigArgs$use.cholesky<- LKrigObject$Mc
 LKrigArgs$wU<- LKrigObject$wU
#
  capture.evaluations <-  rbind( c(lambdaTemp,
                                 a.wghtTemp,
                                 LKrigObject$rho.MLE.FULL,
                                 LKrigObject$sigma.MLE.FULL,
                                 LKrigObject$lnProfileLike.FULL) 
                                )
  if(verbose){
    cat("capture.evaluations first call", fill=TRUE )
    cat( capture.evaluations, fill=TRUE)
  }
  
#####  optimze likelihood over log lambda  and over omega =  log( a.wght -4)/2
    capture.env <- environment()
# last two arguments are specific to this objectinve function     
    result <- try(optim(c(llambda.start, omega.start),
                       LambdaAwghtObjectiveFunction, 
                      lower=c(lowerBoundLogLambda,lowerBoundOmega), 
                      upper=c(upperBoundLogLambda,upperBoundOmega), 
                      method="L-BFGS-B",
                      control=list(fnscale = -1,
                                     ndeps = c(.05,.05),
                                    factr =  1e10),
                      LKrigArgs=LKrigArgs,
                      capture.env= capture.env,
                      verbose=verbose
                      ))
    if(verbose){
      cat("Results from optimize:", fill=TRUE)
      print( result )
    }
    evalSummary <- !(class( result)== "try-error")
    llambda.MLE <- result$par[1]
    lambda.MLE<- exp( llambda.MLE)
    omega.MLE<- result$par[2]
    a.wght.MLE<- 4 + exp( 2 * omega.MLE) 
    LKrigArgs$NtrA <- 20
    
    LKrigObject <- do.call("LKrig", c(LKrigArgs,
                                      list(
                                        lambda = lambda.MLE,
                                        a.wght = a.wght.MLE)
                                     )
    )
  
  ###### end optimze block    
  # save summary results from this set of parameters.
  # Output to be saved     
  out <- rep(NA, 11)
  names( out) <-  c("EffDf", "lnProfLike", "GCV", "sigma.MLE", "rho.MLE", 
                    "lambda.MLE", "a.wght.MLE", "lnLike", "functionEval", 
                    "gradientEval", "totalEval")
    
  out[ 1] <- LKrigObject$trA.est
  out[ 2] <- LKrigObject$lnProfileLike.FULL
  out[ 3] <- LKrigObject$GCV
  out[ 4] <- LKrigObject$sigma.MLE.FULL
  out[ 5] <- LKrigObject$rho.MLE.FULL
  out[ 6] <- lambda.MLE
  out[ 7] <- a.wght.MLE
  out[ 8] <- LKrigObject$lnLike.FULL
  out[ 9] <- result$counts[1]
  out[10] <- result$counts[2]
  out[11] <- nrow(capture.evaluations )
  
  # Name columns  of likelihood eval. 
  dimnames(capture.evaluations)<- list( NULL,
                c("lambda","a.wght", "rho.MLE", "sigma.MLE","lnProfileLike.FULL"))
  return(list(summary = out,
              LKinfo = LKrigObject$LKinfo,
              llambda.start = llambda.start,
              Awght.start=4 + exp( 2 * omega.start),
              lambda.MLE = lambda.MLE,
              a.wght.MLE = a.wght.MLE,
              omega.MLE =omega.MLE,
              llambda.MLE=llambda.MLE,
              lnLike.eval = capture.evaluations,
              call = match.call() )
       )        
}

# Define the objective function 
LambdaAwghtObjectiveFunction<- function(PARS, LKrigArgs, capture.env, verbose=FALSE ) {
  lambdaTemp <- exp( PARS[1] )
  a.wghtTemp <-  omega2Awght( PARS[2], LKrigArgs$LKinfo)
  hold <- do.call("LKrig",          
                  c(LKrigArgs, list(
                    lambda = lambdaTemp,
                    a.wght = a.wghtTemp)
                  )
  )[c("rho.MLE.FULL","sigma.MLE.FULL","lnProfileLike.FULL")] 
  rowForCapture <-c( lambdaTemp,
                     a.wghtTemp,
                     hold$rho.MLE.FULL,
                     hold$sigma.MLE.FULL,
                     hold$lnProfileLike.FULL 
  )
  if( verbose){
    cat(PARS, rowForCapture[c(5,1,2,5)], fill=TRUE )
  }
  lnProfileLike.FULL<- hold$lnProfileLike.FULL 
  temp.eval <- get("capture.evaluations",
                      envir = capture.env )
  assign("capture.evaluations",rbind(temp.eval, rowForCapture),
                      envir = capture.env)
  return(lnProfileLike.FULL)
}

