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

\name{LKrigSAR}
\alias{LKrigSAR}
\alias{LKrigSAR.LKBox}
\alias{LKrigSAR.LKRectangle}
\alias{LKrigSAR.LKInterval}
\alias{LKrigSAR.LKRing}
\alias{LKrigSAR.LKCylinder}
\alias{LKrigSAR.default}
\alias{LKrigSAR.LKSphere}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Method that creates the spatial autoregressive (SAR) matrix. 
	%%  ~~function to do ... ~~
}
\description{Using the information in LKinfo create the SAR matrix 
	for a given level of the multiresolution.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LKrigSAR(object, ...)
  \method{LKrigSAR}{default}(object, ...)
  \method{LKrigSAR}{LKInterval}(object, Level, ...)
  \method{LKrigSAR}{LKRectangle}(object, Level, ...)
  \method{LKrigSAR}{LKBox}(object, Level, ...)
  \method{LKrigSAR}{LKRing}(object, Level, ...) 
  \method{LKrigSAR}{LKCylinder}(object, Level, ...)
   \method{LKrigSAR}{LKSphere}(object, Level, ...)
}
\arguments{
  \item{object}{
  	An LKinfo object. 
}
 \item{Level}{
  	The level of the multiresolution.  
}
  \item{\dots}{
  	Any additional arguments to pass to this method. 
}
}
\details{
	The model for the Gaussian Markov Random field, c, at a given level
	is  
	
	B c  =  e,
	
	where B is the SAR matrix computed by this method, 
	 and e are uncorrelated N(0,1). 
	The precision matrix for this level is  
	
	Q= t(B)\%*\% B 
	
	and so the covariance matrix for c is the inverse of Q:
	
	solve( Q)= solve(B)\%*\% t( solve(B))
	
  
}
\value{
	A matrix in the sparse matrix format, \code{spind}, with dimensions given by the 
	number of lattice points at \code{Level}.
		Because this construction is geometry dependent the default version of
	this method just returns an error message.
}

\author{
	Doug Nychka

}
\seealso{
 \code{\link{LKrig.precision}} \code{\link{spind2full}}
}
\examples{
    x<- cbind( c(0,1))
	LKinfo<- LKrigSetup(x,LKGeometry="LKInterval",
	               nlevel=3, NC=3, a.wght=5, alpha=c(1,.5,.2) )
	B<- LKrigSAR( LKinfo, Level=2)
	B<-spind2full(B)
	image.plot( B)
	
	LKinfo<- LKrigSetup(cbind( c(0,360), c(0,1)) ,LKGeometry="LKRing",
	               nlevel=1, NC=3, a.wght=5, alpha=1)
	B<- LKrigSAR( LKinfo, Level=1)
	B<-spind2full(B)
	image.plot( B)
}
\keyword{spatial}


