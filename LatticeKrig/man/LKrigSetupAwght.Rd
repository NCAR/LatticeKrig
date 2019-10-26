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

\name{LKrigSetupAwght}
\alias{LKrigSetupAwght}
\alias{LKrigSetupAwghtObject}
\alias{LKrigSetupAwght.default}
\alias{LKrigSetupAwght.LKRectangle}

\title{Method to create a.wght component from the \code{LKinfo}
	object.
%%  ~~function to do ... ~~
}
\description{This method takes a vector or more complicated object and
	based on the geometry creates a list with the a.wght information.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LKrigSetupAwght(object, ...)
\method{LKrigSetupAwght}{default}(object, ...)

\method{LKrigSetupAwght}{LKRectangle}(object, ...)

LKrigSetupAwghtObject(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  		The partial or complete \code{LKinfo} object created within \code{LKrigSetup}.
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{
  	Any additional arguments to this method.
%%     ~~Describe \code{\dots} here~~
}
}
\details{
	The simplest function of this method is to convert the 
 	a.wght value into a list that has the length of the number of levels.
 	If only a scalar a.wght value is supplied then the default method
 	just repeats this for each level. 
 	
 	The function \code{LKrigSetupAwghtObject} uses the
 	a.wghtObject component  in the LKinfo object to fill
 	in a.wght parameters for the different
 	levels. This is convenient because the lattice locations are different at each level. The parameters are filled in at level, \code{Level} according to 
 	
 	\preformatted{
 	 latticeLocations<- make.surface.grid( 
 	                          object$latticeInfo$grid[[Level]])
 	 a.wght<- predict( object$a.wghtObject, latticeLocations )
 	}
 	here the predict function is whatever is supplied according to the 
 	class for a.wghtObject. Note that since the returned set of 
 	parameters will be in the format used internally a.wght here will be
 	a list with each component being a matrix. Number of rows are each to
 	the number of lattice points (or basis functions) at that level. 
 	This is easier implement that it may seem and see the examples in \link{nonstationaryModels}.
 	
 	The attribute
 	\code{fastNormalize} (either TRUE or FALSE) is attached to this 
 	list to indicate how the marginal variance of the process should be
 	found.
 	  
\preformatted{ 	 	
  LKinfo<- LKrigSetup( x,LKGeometry="LKInterval", alpha=c( 1,.2,.01),
                   nlevel=3, a.wght=4.5, NC= 3)
  LKrigSetupAwght( LKinfo)
  
[[1]]
[1] 4.5

[[2]]
[1] 4.5

[[3]]
[1] 4.5

attr(,"fastNormalize")
[1] FALSE
}                   
 	Currently the only geometry with fastNormalization being TRUE is for a
 	rectangular domain. 
 	
 	 For the LKRectangle geometry, however, more complicated
	anisotropic and non-stationary a.wght specifications are possible. 
	See \code{\link{LKrig}} for details. Also in the case that the
    fastNormalization is TRUE for rectangles several more attributes are added to the  \code{a.wght} list 	that precompute some matrices of the SAR. 
	
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A list with \code{nlevel} components. The attribute \code{fastNormalize}
	is added to this list. In the case that the geometry is LKRectangle several
	more attributes are added indicating the type of covariance model and
	possibly an eigen decomposition of the SAR matrix exploiting Kronecker
	products.  
}

\author{
	Doug Nychka
%%  ~~who you are~~
}


\seealso{
	\code{\link{LKrigSetup}}, 	\code{\link{LKrigSetupAlpha}},
	\code{\link{LKrigSAR}},  \code{\link{LKrig}}

}
\examples{
  x<- cbind( c(0,1))
  LKinfo<- LKrigSetup( x,LKGeometry="LKInterval", alpha=c( 1,.2,.01),
                   nlevel=3, a.wght=4.5, NC= 3)
  a.wghtList<- LKrigSetupAwght( LKinfo)
  
  x<- cbind( c(0,1), c(0,1))
  LKinfo<- LKrigSetup( x, alpha=c( 1,.2,.01),
                   nlevel=3, a.wght=4.5, NC= 3)
  a.wghtList<- LKrigSetupAwght( LKinfo)
# see   
  names(attributes( a.wghtList))
 
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }


