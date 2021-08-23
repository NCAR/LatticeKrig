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

LKrig.coef <- function(GCholesky, wX, wU, wy, lambda,
                       collapseFixedEffect = FALSE,
                       verbose=FALSE) {
    if (length(lambda) > 1) {
        stop("lambda must be a scalar")
    }
    
    if( !is.null(wU) ){
 #       cat( "HERE1",fill=TRUE)
        A <- forwardsolve(GCholesky, transpose = TRUE, t(wX) %*% wU, 
                              upper.tri = TRUE)
        
        A <- backsolve(GCholesky, A)
        
        A <- t(wU) %*% (wU - wX %*% A)/lambda
#   A is  (U^t M^{-1} U)
        b <- forwardsolve(GCholesky, transpose = TRUE, t(wX) %*% wy, upper.tri = TRUE)
        b <- backsolve(GCholesky, b)
        b <- t(wU) %*% (wy - wX %*% b)/lambda
# b is   (U^t M^{-1} y)
# Save the intermediate matrix   (U^t M^{-1} U) ^{-1}
# this proportional to the GLS covariance matrix of estimated coefficients
# should be small for this to be efficient code 
        hold<- svd( A)
 #       cat( "HERE2",fill=TRUE)
        if( max(hold$d)/min(hold$d) > 1e15){
            print( hold$d )
            stop("large condition number (> 1e15)  in fixed part (X) of model")
        }
 #       cat( "HERE3",fill=TRUE)
        Omega <- solve(A)
 #       cat( "HERE4",fill=TRUE)
# GLS  and also the spatial process estimates for fixed linear part of the model
        d.coef <- Omega %*% b
# combine the different fixed effects estimates across replicates.  
        if(  collapseFixedEffect ){
           d.coefMean<- rowMeans( d.coef)
           dimTemp<- dim ( d.coef)
           d.coef<- matrix( d.coefMean,
                            nrow = dimTemp[1],
                            ncol = dimTemp[2])
        }
        residualFixed<- wy - wU %*% d.coef
   }
   else{
       Omega<- NULL
       d.coef<- NULL
       residualFixed<- wy
   }     
# coefficients of basis functions.
#  cat( "HERE5",fill=TRUE)
    c.coef <- forwardsolve(GCholesky, transpose = TRUE,
                       t(wX) %*% (residualFixed), upper.tri = TRUE)
# Next  is the (very strange) formula from the LKrig article to 
# to evaluate  r^T M^{-1} r  where r = y - U d.coef i.e. r are the residuals from fixed part of the model
# The W from this formula is absorbed into weighting of
# the observations and fixed part of the model. ( e.g. wy is  W^{1/2}y )
# Note that older versions have a mistake 
# where the factor (1/lambda) has been omitted. quad.form is not needed for the 
# coefficients but is used in computing the likelihood.
#   
#  Why not just compute  r^T M^{-1} r directly?  Direct computation requires the 
#  dense and large covariance matrix of the process and so is prohibitive.
#  This indirect form is the results of the magic from the Sherman-Morrison-Woodbury matrix identity.
    
   quad.form<-  (1/lambda) * c( colSums(as.matrix(residualFixed^2))  - 
                                    colSums( as.matrix( c.coef^2) ) )
   #cat( "HERE6",fill=TRUE)
   c.coef <- backsolve(GCholesky, c.coef)
  # cat( "HERE7",fill=TRUE)
   if( verbose){
    	cat("d.coef: ", d.coef, fill=TRUE)
    	cat( fill=TRUE)
    	cat("c.coef: ", c.coef, fill=TRUE)
    }
    return( list(c.coef = c.coef, d.coef = d.coef,
                 Omega = Omega, quad.form=quad.form,
                 collapseFixedEffect= collapseFixedEffect )
            )
}

