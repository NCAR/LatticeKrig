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

\name{setDefaultsLKinfo}
\alias{setDefaultsLKinfo}
\alias{setDefaultsLKinfo.default}
\alias{setDefaultsLKinfo.LKRectangle}
\alias{setDefaultsLKinfo.LKBox}
\alias{setDefaultsLKinfo.LKRing}
\alias{setDefaultsLKinfo.LKCylinder}
\alias{setDefaultsLKinfo.LKInterval}
\alias{setDefaultsLKinfo.LKSphere}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Method for including default information in the LKinfo object.
%%  ~~function to do ... ~~
}
\description{This  method is used to define the various parts of the 
	LatticeKrig model by including specific components in the LKinfo object. 
	The specific instance of the function that is used depends on the
	geometry that has been specified. Typically one would would not need
	to use this
    method for a new geometry. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
setDefaultsLKinfo(object, ...)
\method{setDefaultsLKinfo}{default}(object, ...)
\method{setDefaultsLKinfo}{LKRectangle}(object, ...)
\method{setDefaultsLKinfo}{LKBox}(object, ...)
\method{setDefaultsLKinfo}{LKRing}(object, ...)
\method{setDefaultsLKinfo}{LKCylinder}(object, ...)
\method{setDefaultsLKinfo}{LKInterval}(object, ...)
\method{setDefaultsLKinfo}{LKSphere}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  	An object of class LKinfo.
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{Additional arguments to a specific method. 
%%     ~~Describe \code{\dots} here~~
}
}
\details{
	This method is used to include some default settings or components  in the
	 LKinfo  
	object and is called from \code{LKrigSetup} with the initial (or partial) LKinfo list
	that is formed from the
	information passed by the user to LKrigSetup.
	
	In creating a new geometry this method is not required but can be used to 
	include convenient default values and also some checks on the arguments passed to  LKrigSetup.
	Often there are several arguments for a geometry that make sense to set to simplify the use or
	to avoid bad things .... 
	
	For the LKRing geometry the defaults values follow the logic:
\itemize{	
	\item alpha set to  1.0 if not specified and nlevel is 1 
	
	\item a.wght set to 4.01 if not specified. 
	
     \item the fixed part of the model (fixedFunction) uses the function LKrigPeriodicFixedFunction which insures the first component is periodic.

     \item the transformation matrix, V has the default  diag( c(1,1)). Also in this
	function if V[1,1] is not 1 then an error is given.
	}
	 	For the LKSphere geometry the defaults values follow the logic:
\itemize{	
	\item alpha set to  1.0 if not specified and nlevel is 1 
	
	\item a.wght set to 6.01 if not specified. 
	
     \item the fixed part of the model (fixedFunction) just includes a constant ( m=1)

    \item Default distance function is great circle distance, with sphere radius of 1.0. 
	}
	 

}

\author{
	Doug Nychka
}

\seealso{
	\code{\link{LKrigSetup}}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial }

