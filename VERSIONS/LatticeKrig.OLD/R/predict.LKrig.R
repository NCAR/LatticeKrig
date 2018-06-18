# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2012
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

predict.LKrig <- function(object, xnew = NULL, Znew = NULL, 
    drop.Z = FALSE,return.levels=FALSE, ...) {
newXPassed<- !is.null(xnew)
   if (!newXPassed) {
        xnew <- object$x    
    }
includeZCovariate<- !drop.Z & (object$nZ > 0)
NG <- nrow(xnew)
# setting up Z covariate and some checks
if(includeZCovariate){
    if( newXPassed & is.null(Znew)){
      stop("xnew has been specified, but no new Z covariates ")}
    if (is.null(Znew) ) {
        Znew <- object$Z }
# sanity check on xnew and Znew 
    if( NG!= nrow( Znew)){
         stop(" x (locations) and Z (covariates) have different numbers of
rows")}
     }
# predictions for fixed part of the model
# and can be with or without the additional covariates, Znew.
    T.matrix <- LKrig.fixed.component(xnew,  m = 2, 
        distance.type = object$LKinfo$distance.type)
    if (includeZCovariate) {
       temp1 <- cbind(T.matrix,Znew) %*% object$d.coef
    }
    else {
# in the case of additional covariates only use the coefficients associated
# with the spatial drift 
       temp1 <- T.matrix %*% object$d.coef[object$ind.drift, ]
    }
    PHIg <- LKrig.basis(xnew, object$LKinfo)
# the nonparametric component from the spatial process
# described by the multiresolution basis
  if( !return.levels){
    temp2 <- PHIg %*% object$c.coef
    return(temp1 + temp2)}  
  else{
    nLevels<- object$LKinfo$nlevel
    temp2<- matrix( NA, ncol=nLevels, nrow= nrow( xnew))
    for( level in 1:nLevels){
# indices for each multiresolution level      
       startLevel<- object$LKinfo$offset[level] +1
       endLevel<-  object$LKinfo$offset[level+1]
       indexLevel<- startLevel : endLevel
       print( range( indexLevel))
       temp2[,level]<- PHIg[,indexLevel] %*% object$c.coef[indexLevel]
    }
    temp2[,1] <- temp2[,1] 
    return( cbind(temp1,temp2))
  }
}

