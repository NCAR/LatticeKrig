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

\name{gridList-class}
\Rdversion{1.1}
\docType{class}
\alias{gridList-class}
\alias{gridListInfo}
\alias{gridList}
\title{Class \code{"gridList"}. A  description of a regular and 
	multidimensional grid.}
\description{
	This object is mainly designed to work with methods that take a 
	set of locations organized on a grid. The object is a list where
	there are as many components as dimensions and  each list component is a vector of values being the grid points in that dimension. It is consistent with the older use of the older \code{grid.list} format used in the fields package.
	This form is somewhat redundant because for an equally spaced grid all one needs is the beginning value, spacing and number of points but it makes it simpler to pass the grid information to functions such as image and contour. 
	
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\usage{
gridListInfo(gridList)
}
\arguments{
  \item{gridList}{A gridList object.}
}

\section{Objects from the Class}{

This object is a list where each component is a vector of 
grid points in a particular dimension. For example

\code{ grid<- structure(list( x= seq( -1,1,,20), y= seq( 0,1,,15)), class= "gridList"} 

would create this object for a 2d grid with 20 and 15 points over the 
ranges [-1,1] and [0,1]. 
The component names ( "x" and "y" in this case) are optional. 

The function \code{gridListInfo} extracts some summary information that is used to support the summary function for this class. 
%%  ~~ describe objects here ~~ 
}
\section{Methods}{
  \describe{
    \item{LKrigDistance}{\code{signature(x1 = "matrix", x2 = "gridList", delta = "numeric")}: ... }
	 }
}
\author{
	Doug Nychka
%%  ~~who you are~~
}
\seealso{
	LKrigDistance  and  LKrigDistGrid LKrigLatticeCenters
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("gridList")
# a 3-d grid
grid<- structure(
 list( x= seq( -1,1,,20), y= seq( 0,1,,15) ,oneMore = 1:10) ,
 class= "gridList" )
}
\keyword{classes}


