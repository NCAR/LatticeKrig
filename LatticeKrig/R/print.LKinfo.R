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

print.LKinfo <- function(x, ...) {
#  print( x$call)
    LKinfo <- x
    L <- LKinfo$nlevel
    cat("Classes for this object are: " , class( LKinfo), fill=TRUE)
    
    cat("The second class usually will indicate the geometry
     e.g.  2-d rectangle is  LKRectangle", fill=TRUE)
    cat(" ", fill = TRUE) 
    cat("Some details on spatial autoregression flags:",     fill=TRUE)
    cat( "stationary: ",  attr( LKinfo$a.wght,"stationary"),  fill=TRUE)
    cat("first order (by level): ",  attr( LKinfo$a.wght,"first.order"), fill=TRUE)
    cat("isotropic: "  ,   attr( LKinfo$a.wght,"isotropic"),   fill=TRUE)
    cat(" ", fill = TRUE)
    if(LKinfo$dense){
      cat("Hey! The dense flag is TRUE so computations will
          not be done using sparse matrices.", fill=TRUE)
    }
    cat("Ranges of locations in raw scale:", fill=TRUE)
    print(  LKinfo$latticeInfo$rangeLocations)
    if( !is.null(LKinfo$basisInfo$V)){
    	cat("(inverse) linear transformation for lattice nodes:",fill=TRUE)
    	print(LKinfo$basisInfo$V )
    	cat("transformed ranges:",fill=TRUE)
    	print( LKinfo$latticeInfo$grid.info$range)
    }
    cat(" ", fill = TRUE)
    cat("Logical (collapseFixedEffect) if fixed effects will be pooled:" ,
        LKinfo$collapseFixedEffect, fill = TRUE)
    cat(" ", fill = TRUE)
    cat("Number of levels:", L, fill = TRUE)
    cat("delta scalings:", x$latticeInfo$delta, fill = TRUE)
    cat("with an overlap parameter of ", LKinfo$basisInfo$overlap, fill=TRUE)
    cat("alpha: ", unlist(x$alpha), fill = TRUE)
    if (!is.null(x$nu)) {
        cat("based on smoothness nu = ", x$nu, fill = TRUE)
    }
    if(!is.null( x$alphaObject)){ 
      cat("alpha specified at each level by objects: ",fill=TRUE)
      cat("Level ", "Class", fill=TRUE)
      for( level in 1:L){
               cat( level, class(x$alphaObject[[1]]), fill=TRUE )
      }
    }
    temp<- unlist(x$a.wght)
    
    cat(" ", fill = TRUE)
    if( length( temp)<= x$nlevel * 9){
    cat("a.wght: ", temp, fill = TRUE)
    }
    else{
      cat("dim(A.wght[[k]]): ",
          fill=TRUE )
      for ( k in 1: x$nlevel){
      cat(" Level", k,
          dim(x$a.wght[[k]]), fill=TRUE )
      }
    }
    if(!is.null( x$a.wghtObject)){ 
      cat("a.wght specified at each level by objects: ",fill=TRUE)
      cat("Level ", "Class","subClass", fill=TRUE)
      for( level in 1:L){
        cat( level, class( x$a.wghtObject), class(x$a.wghtObject[[1]]), fill=TRUE )
      }
    }
       cat(" ", fill = TRUE)
  # Details on basis functions at each level
      bType <- LKinfo$basisInfo$BasisType
      cat( "Basis  type:",
           LKinfo$basisInfo$BasisType, 
           "using ",
           LKinfo$basisInfo$BasisFunction,
           " and", LKinfo$distance.type, " distance.", fill=TRUE)
               if( LKinfo$normalize){
  cat("Basis functions will be normalized", fill=TRUE)
        }
  cat(" ", fill = TRUE)      
  cat("Total number of basis functions ",  LKinfo$latticeInfo$m, fill=TRUE)  
      temp<- cbind( 1:L, LKinfo$latticeInfo$mLevel)
      dimnames( temp)<- list( rep( "", L), c("Level" ,   "Basis size"))
     if( !is.null(LKinfo$latticeInfo$mx)){
    	 temp<- cbind( temp, LKinfo$latticeInfo$mx)
 	 }
      print( temp)      	
        cat(" ", fill = TRUE)  
 cat("Lambda value: ", LKinfo$lambda, fill=TRUE)         
}



