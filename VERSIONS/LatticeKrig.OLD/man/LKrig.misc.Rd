% # LatticeKrig  is a package for analysis of spatial data written for
% # the R software environment .
% # Copyright (C) 2012
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

% # Modified DH Dec 2013
% Earlier modifications DH 2/03/2012

\name{LKrig Miscellaneous Matrix Functions}
\alias{LKrig.rowshift.periodic}
\alias{ LKrig.shift.matrix}
\alias{LKrig.rowshift}
\alias{which.max.matrix}
\alias{which.max.image}
\title{
Miscellaneous internal functions for LatticeKrig package. 
}
\description{
Some utility functions used internally by higher level \code{LKrig}
functions. Currently these are simple functions that perform shifts of
a matrix.
}
\usage{
LKrig.shift.matrix( A, shift.row=0, shift.col=0, periodic=c(FALSE, FALSE))
LKrig.rowshift.periodic( A, shift.row)
LKrig.rowshift( A,shift.row,shift.col)
which.max.matrix(z)
which.max.image(obj)
}

\arguments{

\item{A}{A matrix.}


\item{obj}{An image list with the usual components x, y, and z.}

\item{periodic}{ A vector of two logicals pertaining to rows and
columns. TRUE indicates an index where the shift will be periodic --
entries shifted beyond the dimensions will be wrapped to the other
side.}

\item{shift.row}{An integer that specifies the number of positions
 that the rows of the matrix are shifted.}

\item{shift.col}{An integer that specifies the number of positions
 that the columns of the matrix are shifted.}

\item{z}{A matrix.}

} \details{

These functions are used to create the nearest neighbor indices for
the precision matrices.

}
\value{

A matrix of shifted values. Entries that are not defined due to
the shift are set to NA.  A column shift is done by a combination of
transpose operations and a row shift.

\preformatted{
A<- matrix( 1:12,3,4)
A
     [,1] [,2] [,3] [,4]
[1,]    1    4    7   10
[2,]    2    5    8   11
[3,]    3    6    9   12

#shift of 2 for rows:
 LKrig.rowshift( A, 2)
    [,1] [,2] [,3] [,4]
[1,]   NA   NA   NA   NA
[2,]   NA   NA   NA   NA
[3,]    1    4    7   10

#periodic case
LKrig.rowshift.periodic( A, 2)
     [,1] [,2] [,3] [,4]
[1,]    2    5    8   11
[2,]    3    6    9   12
[3,]    1    4    7   10
}

}

\author{
Doug Nychka
}

\seealso{
LKrig.MRF.precision
}

\keyword{spatial}
