
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

\name{LKrigSetupLattice}
\alias{LKrigSetupLattice}
\alias{LKrigSetupLattice.default}
\alias{LKrigSetupLattice.LKBox}
\alias{LKrigSetupLattice.LKInterval}
\alias{LKrigSetupLattice.LKRectangle}
\alias{LKrigSetupLattice.LKRing}
\alias{LKrigSetupLattice.LKCylinder}
\alias{LKrigSetupLattice.LKSphere}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates the lattice information for a specific geometry.
%%  ~~function to do ... ~~
}
\description{Given a specific geometry and the initial information supplied in an 
	LKinfo list create the information that is needed to define the lattice for a
	given model. This function is required for any new geometry added to LatticeKrig.  
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LKrigSetupLattice(object, ...)
\method{LKrigSetupLattice}{default}(object, ...)

\method{LKrigSetupLattice}{LKBox}(object, verbose,...)
 
\method{LKrigSetupLattice}{LKRectangle}(object, verbose,...)

\method{LKrigSetupLattice}{LKInterval}(object, verbose,...)

 \method{LKrigSetupLattice}{LKRing}(object, verbose,...)
 
\method{LKrigSetupLattice}{LKCylinder}(object, verbose,...)

 \method{LKrigSetupLattice}{LKSphere}(object, x = NULL, verbose,...)
    
}

\arguments{
%\item{grid.info}{}
  \item{object}{ A list that is an LKinfo object.}


\item{verbose}{If TRUE print out intermediate information for debugging.
}

\item{x}{Locations of the observations that define spatial domain. For LKSphere locations
are in lon/lat degrees.}

\item{\dots}{Any additional arguments. }
  
}
\details{
	This method takes the LKInfo object and the other arguments and computes the 
	lattice information needed for a specific geometry. The functions are
	 called from within LKrigSetup and the results are added as a component to\code{latticeInfo}
	as part of the LKinfo object.  The  way 
	to design  what should be in \code{latticeInfo} is to keep in mind that
	creating the spatial AR matrix (LKrigSAR) and  determining the multi-resolution
	 lattice points (LKrigLatticeCenters) use the information in the LKinfo
	  object. 
 
 Because the lattice must depend on the geometry the default method just
 prints out an error message.
 
   
  \strong{NC  and  NC.buffer} 
  for Cartesian geometries (e.g. \code{LKInterval, LKRectangle, LKBox, 
  LKRing}) are parameters to specify to number of the
  	 grid points in the largest dimension of the spatial domain and for the 
  	 coarsest lattice. These should already be a component in the LKinfo object,
  	 object specified in the usage. 
  	 
	 For \code{LKInterval} this is just of lattice points at the coarsest level.
	 For LKRectangle and LKBox if the spatial domain has the same size in all 
	 dimensions then the number of lattice points within the spatial domain are
	  NC^2 and NC^3 respectively. Note that the total number of lattice points
	   at a given level is also effected the size of NC.buffer (see below). 

  The number of lattice points to add to the margins beyond the
  	 spatial domain is controlled by \code{NC.buffer}. Thus in the largest spatial dimension there are a total of 
		NC + NC.buffer*2 grid points. This is the number of lattice points for 
		LKInterval.
}
\value{
	For the \code{LKInterval, LKRectangle, LKBox} geometries.
	A list with  required components:
	\describe{
		\item{m}{The total number of lattice points. i.e. the total number
			      of basis functions.}
			      \item{mx}{A matrix giving the number of lattice points in 
			      	each coordinate (columns) and at each level (rows) of the 
			      	multi-resolution.}
		\item{offset}{When the lattice points are unrolled as a single array,
			 a vector of indexes giving the start of each 
			 lattice level in the coefficient vector. }
		\item{delta}{A vector of the  lattice point spacings for each level.}
		\item{rangeLocations}{Limits of spatial domain.}
		\item{mLevel}{A vector giving number of lattice points at each level.}
	}
	
	In addition, the methods for Cartesian and cartesian-like  spatial domains (
	LKInterval, LKRectangle, LKBox, LKRing, LKCylinder, LKSphere) include the additional components:
	\describe{	
		\item{mx}{A matrix giving number of grid points at each level and for each dimension.}
		\item{mLevelDomain}{Same as mLevel but restricted to lattice points
			 within the spatial domain.}
		\item{mxDomain}{Same as mx but restricted to the points within the
			 spatial domain.}
		
		\item{NC}{Passed value.}
		\item{NC.buffer}{Passed value.}
		\item{grid}{A list where each component is also list giving the grid points of the lattice in each coordinate. }
    }
    
  LKSphere also returns the additional component \code{grid3d} -- the direction Cosines
  of the grid points. 
    
    
      

}

\author{
	Doug  Nychka
%%  ~~who you are~~
}

\seealso{
	\code{\link{LKGeometry}}
    \code{\link{LKrigSetup}},
	\code{\link{LKrigSAR}},
	\code{\link{LKrigLatticeCenters}}
}


\keyword{ spatial }


