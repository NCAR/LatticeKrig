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

\name{LKrigSetupAlpha}
\alias{LKrigSetupAlpha}
\alias{LKrigSetupAlpha.default}
\alias{LKrigSetupAlpha.LKInterval}
\alias{LKrigSetupAlpha.LKRectangle}
\alias{LKrigSetupAlpha.LKBox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates the alpha parameter list in LatticeKrig covariance.
%%  ~~function to do ... ~~
}
\description{This function is called by \code{LKrigSetup} and
	creates the list for the alpha parameters based on the information 
	from the call to \code{LKrigSetup} and additional information and 
	conditions related to the geometry. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LKrigSetupAlpha(object, ...)
\method{LKrigSetupAlpha}{default}(object, ...)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  	The partial LKinfo object created within \code{LKrigSetup}
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{
  	Any additional arguments to this method
%%     ~~Describe \code{\dots} here~~
}
}
\details{
	The main function of this method is to fill in the sequence of alpha values
	for a parametric model and convert those values to a list instead of a  vector. In the case that the scalar \code{nu} is supplied it 	is used to create the list according to:
\preformatted{	
	alpha <- exp(-2 * (1:nlevel) * nu)
    alpha <- alpha/sum(alpha)
    as.list( alpha)
}    
    
	%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with \code{nlevel} components each  representing the alpha values
at that level. In the simplest case a vector of alpha values is converted to 
a list.  
\preformatted{
  LKinfo<- LKrigSetup( x, alpha=c( 1,.2,.01),
                   nlevel=3, a.wght=4.5, NC= 3)
  LKrigSetupAlpha( LKinfo)
[[1]]
[1] 1
[[2]]
[1] 0.2
[[3]]
[1] 0.01
}

}
\author{
	Doug Nychka

}
\seealso{
	\code{\link{LKrigSetup}}, 	\code{\link{LKrigSetupAwght}},
	\code{\link{LKrigSAR}},  \code{\link{LKrig}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# an x that is just the limits of the domain	
  x<- cbind( c(0,1), c(0,1))
  
  LKinfo<- LKrigSetup( x, alpha=c( 1,.2,.01),
                   nlevel=3, a.wght=4.5, NC= 3)
  alphaList<- LKrigSetupAlpha( LKinfo)

  LKinfo<- LKrigSetup( x, nu=1, nlevel=4, a.wght=4.5, NC= 4)
  alphaList<- LKrigSetupAlpha( LKinfo)
}

\keyword{spatial}
% __ONLY ONE__ keyword per line

