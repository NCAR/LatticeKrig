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

LKrigMarginalVariance<- function(x1, LKinfo, verbose = FALSE)
  {
    nlevel        <- LKinfo$nlevel
    delta         <- LKinfo$latticeInfo$delta
    overlap       <- LKinfo$basisInfo$overlap
    normalize     <- LKinfo$normalize
    distance.type <- LKinfo$distance.type
#    fast          <-  attr( LKinfo$a.wght,"fastNormalize")
# do not use fast normalize just ot keep code simple
    V <- LKinfo$basisInfo$V
# coerce x1 to a matrix    
    x1<- as.matrix( x1)
   
    # transform locations if necessary (lattice centers already in 
    # transformed scale) 
    if( !is.null( V[1]) ){     
      x1<- x1 %*% t(solve(V))     
    }
    if( verbose){
          cat("LKrig.basis: Dim x1 ",  dim( x1), fill=TRUE)
    }
    marginalVar<- matrix( NA, ncol=nlevel,  nrow=nrow( x1))
    for (l in 1:nlevel) {
        # Loop over levels and evaluate basis functions in that level.
        basis.delta <- delta[l] * overlap
        # 
        # There are two choices for the type of basis functions
        # 
        centers<- LKrigLatticeCenters( LKinfo,Level=l )
        if(LKinfo$basisInfo$BasisType=="Radial" ){ 
        PHItemp <- Radial.basis(  x1, centers, basis.delta,
                                max.points = LKinfo$basisInfo$max.points,
                             mean.neighbor = LKinfo$basisInfo$mean.neighbor, 
                       BasisFunction = get(LKinfo$basisInfo$BasisFunction),
                             distance.type = LKinfo$distance.type,
                                   verbose = verbose)
                             
                             }
        if(LKinfo$basisInfo$BasisType=="Tensor" ){  
                
        PHItemp <- Tensor.basis(  x1, centers, basis.delta,
                                max.points = LKinfo$basisInfo$max.points,
                             mean.neighbor = LKinfo$basisInfo$mean.neighbor, 
                       BasisFunction = get(LKinfo$basisInfo$BasisFunction),
                             distance.type = LKinfo$distance.type)
                             
                             }      	                            

        	# the default choice should work for all models	
               marginalVar[,l]<-  LKrigNormalizeBasis( LKinfo,  Level=l,  PHI=PHItemp)  
    }
    return(marginalVar)
}


