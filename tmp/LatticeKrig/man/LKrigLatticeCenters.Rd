% # LatticeKrig  is a package for analysis of spatial data written for
% # the R software environment .
% # Copyright (C) 2016
% # University Corporation for Atmospheric Research (UCAR)
% # Contact: Douglas Nychka, nychka@ucar.edu,
% # National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
% #
% # This program is free software; you can redistribute it and/or modify
% # it under the terms of the GNU General Public License as published by
% # the Free Software Foundation; either version 2 of the License, or
% # (at your option) any later version.
% # This program is distributed in the hope that it will be useful,
% # but WITHOUT ANY WARRANTY; without even the implied warranty of
% # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% # GNU General Public License for more details.

\name{LKrigLatticeCenters}
\alias{LKrigLatticeCenters}
\alias{LKrigLatticeCenters.default} 
\alias{LKrigLatticeScales}
\alias{LKrigLatticeScales.default}
\alias{LKrigLatticeCenters.LKBox}        
\alias{LKrigLatticeCenters.LKInterval}     
\alias{LKrigLatticeCenters.LKRectangle}
\alias{LKrigLatticeCenters.LKRing}
\alias{LKrigLatticeCenters.LKCylinder}
\alias{LKrigLatticeCenters.LKSphere}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Methods to report the locations  or scales associated with the lattice points.
%%  ~~function to do ... ~~
}
\description{These method takes the lattice information for a particular geometry 
	from an LKinfo object and finds the locations or scales at each lattice points.
	 These locations are the 
	"nodes" or centers of the basis functions. The "scales" scales that distance function when the 
	basis functions are evaluated and combine the spacing of lattice and the specified overlap. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LKrigLatticeCenters(object, ...)
\method{LKrigLatticeCenters}{default}(object, ...)
\method{LKrigLatticeCenters}{LKInterval}(object, Level, ...)
\method{LKrigLatticeCenters}{LKRectangle}(object, Level, ...)
\method{LKrigLatticeCenters}{LKBox}(object, Level, ...)
\method{LKrigLatticeCenters}{LKCylinder}(object, Level = 1, physicalCoordinates = FALSE, ...)
\method{LKrigLatticeCenters}{LKRing}(object, Level = 1, 
                    physicalCoordinates = FALSE, ...)
\method{LKrigLatticeCenters}{LKSphere}(object, Level, ...)
\method{LKrigLatticeScales}{default}(object, ...)
LKrigLatticeScales(object, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  	An LKinfo object. 
}
\item{Level}{
The multi-resolution level.
}
\item{physicalCoordinates}{If TRUE the centers are returned in the 
untransformed scale. See the explanation of the 
\code{V} matrix in LKrigSetup. This is useful to relate the lattice points to 
other physical components of the problem. 

For example with the LKRing geometry 
representing the equatorial slice of the solar atmosphere one observes a line of 
sight integral through the domain. This integral is obvious found with respect to the physical coordinates and not the lattice points. 
}

  \item{\dots}{
  	Any additional arguments for this method.

}
}
\details{
	This method is of course geometry dependent and the default version just
	gives an error warning that a version based on the geometry is required. 
	Typically generating these lattice points from the information in LKinfo
	should be easy as the grid points are already determined. 
	
	The scales reported are in the simplest form  delta*overlap where delta is a 
	vector of the lattice spacings and overlap (default is 2.5) is the amount of 
	overlap between basis functions. 
	
	See the source for the function \code{\link{LKrig.basis}} for how each of these
	is used to evaluate the basis functions. 
	
	
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\strong{Centers} A matrix where the rows index the points and columns index dimension. 
In the case of the LKRectangle geometry attribute is added to indicate the
grid points used to generate the lattice.  For LKSphere the centers are in lon/lat degrees.
( Use \code{directionCosines} to convert to 3-d coordinates from lon/lat.)

\strong{Scales} The default method returns the vector  \code{delta*offset} with length being the
number of multi-resolution levels. 
}

\author{
	Doug Nychka
%%  ~~who you are~~
}
\seealso{
	    \code{\link{LKrig.basis}}
		\code{\link{LKrigSetup}}, 	\code{\link{LKrigSetupAwght}},
	\code{\link{LKrigSAR}},  \code{\link{LKrig}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  x<- cbind( c(-1,2), c(-1,2))
  LKinfo<- LKrigSetup( x, alpha=c( 1,.2,.01),
                   nlevel=3, a.wght=4.5, NC= 10)
# lattice centers for the second level   
# not points added for buffer outside of spatial domain                
   look<- LKrigLatticeCenters(LKinfo, Level=2) 
# convert grid format (gridList)  to just locations
   look<- make.surface.grid( look)
   plot( look,  cex=.5)
   rect( -1,-1,2,2, border="red4")
   
    x<- cbind( c(0, 360), c( 1,3)) 
    LKinfo<- LKrigSetup( x, LKGeometry="LKRing",
                   nlevel=1, a.wght=4.5, NC= 10, V= diag(c( 1,.01) ) )
                   
    polar2xy<- function(x){
	x[,2]*cbind( cos(pi*x[,1]/180), sin(pi*x[,1]/180))}
	        
    look1<- LKrigLatticeCenters( LKinfo, Level=1)               
    look2<- LKrigLatticeCenters( LKinfo, Level=1, physicalCoordinates=TRUE )
    look3<- polar2xy( look2$Locations )
# Basis scales:    
      LKrigLatticeScales( LKinfo)
    
    set.panel(3,1)
    plot( make.surface.grid( look1))
    plot( look2$Locations)
    plot( look3)

                 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }

