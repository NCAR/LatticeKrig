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

\name{LKrigDistance-methods}
\docType{methods}
\alias{LKrigDistance-methods}
\alias{LKrigDistance}
\alias{LKrigDistance,matrix,gridList,numeric-method}
\alias{LKrigDistance,matrix,matrix,numeric-method}
\title{ Distance function methods for \code{LKrigDistance} in Package \pkg{LatticeKrig} ~~}
\description{Distance function for distances within a 
	limited radius and creating a sparse matrix. Supports either
	coordinates or locations specified by a grid. This method allows for 
	finding pairwise distances among locations where one set can 
 ~~ Methods for function \code{LKrigDistance} in package \pkg{LatticeKrig} ~~
}
\section{Methods}{
\describe{	

\item{\code{signature(x1 = "matrix", x2 = "gridList", delta = "numeric")}}{
Finds pairwise distances within the radius delta between all locations
\code{x1} and all grid points. Will call either \code{LKDistGrid} or \code{LKDistGridComponents}. The function is particularly efficient
when at least one set of locations is a regular grid and this is the 
main reason for this "overloading". 	
%%  ~~describe this method here~~
}

\item{\code{signature(x1 = "matrix", x2 = "matrix", delta = "numeric")}}{
Finds pairwise distances within the radius delta between all locations
\code{x1} and \code{x2}. Will call either \code{LKdist} or \code{LKDistComponents}.
%%  ~~describe this method here~~
}
}}
\keyword{methods}
\keyword{spatial}

