# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2

LKrigFindLambda <- function(x, y, ...,  LKinfo,
                            use.cholesky=NULL, 
                 lambda.profile=TRUE,
            lowerBoundLogLambda=-16,
                            tol=.005,
                        verbose=FALSE) {
# parts of the LKrig call that will be fixed.  (except updates to LKinfo)                           	
    LKrigArgs <- c(list(x = x, y = y), list( ...),
                   list( LKinfo=LKinfo, NtrA=ifelse( lambda.profile, 0, 20) ))
    if( verbose){
    	cat( "LKrigFindLambda: Complete set of LKrigArgs:", names(LKrigArgs ), fill=TRUE)
    }               
   
    capture.evaluations <- matrix(NA, ncol = 4, nrow = 1, 
                dimnames = list(NULL, c("lambda", "rho.MLE",
                                         "sigma.MLE", "lnProfileLike.FULL")))
    optim.counts<-  NA
    llambda.start<- log( LKrigArgs$LKinfo$lambda )
    if( is.na(llambda.start)){ llambda.start<- -1 }
    if(verbose){
    	cat("LKrigFindLambda: llambda.start:",  llambda.start, fill=TRUE)
    }
#    
# first fit to get cholesky symbolic decomposition  and wX and wU matrices
# Note that if use.cholesky != NULL then the symbolic decomposition information is
# used from the passed object.

    LKrigObject <- do.call("LKrig", c(
                             LKrigArgs,
                             list( 
                                   use.cholesky = use.cholesky, 
                                return.cholesky = TRUE,
                                 return.wXandwU = TRUE,
                                         lambda = exp( llambda.start),
                                        verbose = verbose)))
                                   
    capture.evaluations[1,] <-  c( LKrigObject$LKinfo$lambda,
                                   LKrigObject$rho.MLE.FULL,
                                   LKrigObject$sigma.MLE.FULL,
                                   LKrigObject$lnProfileLike.FULL)                                                   
    llambda.opt<- llambda.start
    Mc.save<- LKrigObject$Mc
    wX.save<- LKrigObject$wX
    wU.save<- LKrigObject$wU
#    print( dim( wX.save))
#    
##### in most cases now optimze likelihood over log lambda    
   if( lambda.profile){
#  temporary function used in optimizer
#   
# objective function    
    temp.fn <- function(x) {
#    cat("temp.fn: lambda ", exp(x), fill=TRUE)
         lambdaTemp<- exp(x)                       
         hold <- do.call("LKrig",          
                     c(LKrigArgs, list(
                            use.cholesky = Mc.save,
                                      wX = wX.save,
                                      wU = wU.save,
                                 verbose = FALSE,
                                  lambda = lambdaTemp 
                                  )
                       )
                       )[c( "rho.MLE.FULL", "sigma.MLE.FULL", 
                                   "lnProfileLike.FULL") ]               
            lnProfileLike.FULL<- hold$lnProfileLike.FULL 
            rowForCapture<- c(lambdaTemp,
                              hold$rho.MLE.FULL,
                              hold$sigma.MLE.FULL,
                              hold$lnProfileLike.FULL)
            temp.eval <- get("capture.evaluations")
            assign("capture.evaluations",rbind(temp.eval, rowForCapture),
             envir = capture.env)
                             return(lnProfileLike.FULL)
            }
#    
# the try wrapper captures case when optim fails.   
            capture.env <- environment()
            look<- try(optimize( temp.fn, interval = llambda.start+c(-8,5),
                                   maximum= TRUE, tol=tol))
            if(verbose){
            	cat("Results from optimize:", fill=TRUE)
            	print( look)
            }
            evalSummary <- !(class( look)== "try-error")
            llambda.MLE <- look$maximum  
            lambda.MLE =  exp(llambda.MLE)
            optim.counts <- nrow( capture.evaluations) + 2
            LKrigArgs$NtrA <- 20
# surgery on the args list to update lambda with the MLE.
#            LKrigArgs$LKinfo<- LKinfoUpdate(LKrigArgs$LKinfo,lambda = lambda.MLE  )
            LKrigObject <- do.call("LKrig", c(LKrigArgs,
                                    list(use.cholesky = Mc.save,
                                                   wX = wX.save,
                                                   wU = wU.save,
                                               lambda = lambda.MLE)
                                              )
                                   )
   }
  
###### end optimze block    
# save summary results from this set of parameters.
  if( lambda.profile ){  
   out <-  c( LKrigObject$trA.est,
              LKrigObject$lnProfileLike.FULL,
              LKrigObject$GCV,
              LKrigObject$sigma.MLE.FULL,
              LKrigObject$rho.MLE.FULL,
              lambda.MLE,
              llambda.MLE,
              LKrigObject$lnLike.FULL,
              optim.counts,
              NA)
   
   names( out) <-  c("EffDf", "lnProfLike", "GCV", 
                     "sigma.MLE", "rho.MLE", "lambda.MLE",
                     "llambda.MLE", "lnLike",
                     "counts value", "grad")
   }
   else{
# simple clean up when optimization is not done    
     out <-  c( LKrigObject$trA.est,
                LKrigObject$lnProfileLike.FULL,
                LKrigObject$GCV,
                LKrigObject$sigma.MLE.FULL,
                LKrigObject$rho.MLE.FULL,
                exp( llambda.start),
                llambda.start,
                LKrigObject$lnLike.FULL,
                NA,
                NA)
     
      lambda.MLE<- exp( llambda.start)
      names( out) <-  c("EffDf", "lnProfLike", "GCV", 
                       "sigma.MLE", "rho.MLE", "lambda",
                       "llambda", "lnLike",
                       "counts value", "grad")
   }
   
   return(list(summary = out,
                LKinfo = LKinfo,
         llambda.start = llambda.start,
            lambda.MLE = lambda.MLE,
               lnLike.eval = capture.evaluations
               )
# call omitted because it can substitute actual
# values and get large
#                  call = match.call()
#                  )
          )        

}

