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

\name{LKrigNormalizeBasis}
\alias{LKrigNormalizeBasis}
\alias{LKrigNormalizeBasisFast}                            
\alias{LKrigNormalizeBasisFast.default}
\alias{LKrigNormalizeBasisFast.LKRectangle}
\alias{LKRectangleSetupNormalization}

\title{Methods and functions to support normalizing to 
	marginal variance of one. 

}
\description{The basis functions can be normalized so that the 
	marginal variance of the process at each level and at all locations is
	one. A generic function \code{LKrigNormalizeBasis} computes this for
	any LatticeKrig model. However, in special cases the normalization can 
	be accelerated taking advantage of the geometry and the model for the
	lattice weights. This alternative is the method \code{LKrigNormalizeBasisFast}.

}
\usage{
LKrigNormalizeBasis( LKinfo, Level, PHI, ...)
LKrigNormalizeBasisFast(LKinfo, ...)
\method{LKrigNormalizeBasisFast}{default}(LKinfo, ...)
\method{LKrigNormalizeBasisFast}{LKRectangle}(LKinfo, Level, x,  ...)
LKRectangleSetupNormalization(mx, a.wght)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a.wght}{A.wght parameters.}
  \item{LKinfo}{ An LKinfo object.  NOTE:  Here  \code{choleskyMemory}, a Spam memory argument, can be a component of LKinfo and is subsequently passed
  	 through to the	(spam) sparse cholesky decomposition}
  \item{Level}{The multi-resolution level.}
  \item{mx}{Matrix of lattice sizes.}
  \item{PHI}{Unnormalized basis functions evaluated at the locations 
  	to find the normalization weights.}
  \item{x}{Locations to find normalization weights.}
  \item{\dots}{Additional arguments for method. }
}
\details{
Normalization to unit variance is useful for reducing the artifacts of
the lattice points and creates a model that is closer to being stationary. 
The computation must be done for every point evaluated with the basis functions 
The basic calculation is 
\preformatted{
	tempB <- LKrigSAR(LKinfo, Level = Level)
    tempB <- LKrig.spind2spam(tempB)
    wght <- LKrig.quadraticform(t(tempB) %*% tempB, PHI = PHI, 
        choleskyMemory = choleskyMemory)
}
This generic method uses the low level utility \code{LKrig.quadraticform} that evaluates  diag( t(PHI) solve( Q.l) PHI ) where PHI are the basis functions evaluated at the locations and Q.l is the precision matrix for the lattice coefficients at level l. 
 

For constant a.wght and for the rectangular geometry one can use an eigen decomposition of the Kronecker product of the SAR matrix. 
}
\value{
A vector of weights. The basis are normalized by dividing by the square root
of the weights.  
}
\author{
Doug Nychka
}

\keyword{ spatial }

