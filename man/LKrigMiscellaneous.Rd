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

% # Modified DH Dec 2013
% Earlier modifications DH 2/03/2012

\name{LKrig Miscellaneous Matrix Functions}
\alias{LKrig.rowshift.periodic}
\alias{LKrig.shift.matrix}
\alias{LKArrayShift}
\alias{LKrig.rowshift}
\alias{which.max.matrix}
\alias{which.max.image}
\alias{expandMList}
\alias{expandMatrix0}
\alias{expandMatrix}
\alias{repMatrix}
\alias{convertIndexPeriodic}
\alias{grid2Index}
\alias{convertIndexArray}

\title{
Miscellaneous internal functions for LatticeKrig package. 
}
\description{
Some utility functions used internally by higher level \code{LKrig}
functions. Currently these are simple functions that perform shifts of
a matrix and operations on indices for multidimensional arrays.
}
\usage{
LKrig.shift.matrix( A, shift.row=0, shift.col=0, periodic=c(FALSE, FALSE))
LKrig.rowshift.periodic( A, shift.row)
LKrig.rowshift( A, shift.row, shift.col)
LKArrayShift(A, shift, periodic = FALSE)

which.max.matrix(z)
which.max.image(obj)
expandMatrix0( A, B)
expandMatrix( ...)
expandMList( Mlist, byrow=TRUE)


convertIndexPeriodic(I, nGrid, nPad = NULL)
convertIndexArray(I, nGrid)
grid2Index(I, grid)
}

\arguments{

\item{A}{A matrix.}

\item{byrow}{If TRUE matrices will be repeated row by row. If FALSE
this will be done column by column.}

\item{B}{Another matrix.}
\item{grid}{A vector giving the size of each dimension of array.}

\item{I}{A matrix of multidimensional indices where each row
	identifies an element of the array. e.g. If I[1,] == c( 3,4,2)
	this refers to  A[3,4,2] }

\item{Mlist}{A list where each component is a matrix.}

\item{nGrid}{ A array giving number of elements in each dimension.}

\item{nPad}{An array with the number of padding indices in each dimension.}

\item{obj}{An image list with the usual components x, y, and z.}

\item{periodic}{  A vector of  logicals columns.TRUE indicates an index where the shift will be periodic --
entries shifted beyond the dimensions will be wrapped to the other
side. e.g. for the matrix version c(FALSE, TRUE) indicates the 
columns of the matrix will be shifted in a periodic way but not the 
rows. 
	}
\item{shift}{A integer vector that gives the shifts for each dimension of the
	array.}
\item{shift.row}{An integer that specifies the number of positions
 that the rows of the matrix are shifted.}

\item{shift.col}{An integer that specifies the number of positions
 that the columns of the matrix are shifted.}

\item{z}{A matrix.}

\item{\dots}{Matrices to be expanded.}
} \details{

\strong{Shift related:}
These functions are used to create the nearest neighbor indices for
the precision matrices.

\strong{Expand related:}
These functions are useful for creating a sets of covariance parameters that 
follow a factorial pattern. For example repeating the rows of the "alpha" parameters
as the "a.wght" parameters are varied. \code{expandMList} is particularly useful for creating
a factorial design of parameters to pass to LKrig.MLE for searching the likelihood.

\strong{Index related:}
The function \code{convertIndexPeriodic} converts a  single index for a multidimensional array, into an index that reflects wrapping into a smaller grid to reflect the padding. 
This is used in the \code{LKDistGrid} function to handle distances when the grid has periodic dimensions. The other two functions are used for checking and debugging of going back and
forth between the multidimensional and single indexes.  

}
\value{

\strong{Shift:}
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

\strong{Expand:}
\code{ExpandMList} Returns a list of matrices where the original matrices are repeated so that
combinations of rows are represented. The example below illustrates.  \code{byrow=FALSE} does the repetition
by columns instead of rows.

\preformatted{
> A
     [,1] [,2]
[1,]    1    3
[2,]    2    4
> B
     [,1]
[1,]   11
[2,]   12
[3,]   13
> C
[1,] 100
> expandMList( list( A=A, B=B, C=C))
$A
     [,1] [,2]
[1,]    1    3
[2,]    2    4
[3,]    1    3
[4,]    2    4
[5,]    1    3
[6,]    2    4

$B
     [,1]
[1,]   11
[2,]   11
[3,]   12
[4,]   12
[5,]   13
[6,]   13

$C
     [,1]
[1,]  100
[2,]  100
[3,]  100
[4,]  100
[5,]  100
[6,]  100
}
}

\author{
Doug Nychka
}
\examples{
	A<- array( 1:90, c( 4,5,3))
	LKArrayShift( A, c( -1,-1,0))	
	
# welcome to the world of unrolling multiarray indices
A<- array( 1:60, c( 4,3,5))	
I<- rbind( c(1,2,1), c( 3,2,5))
look<- grid2Index( I, c( 4,3,5) )
# A has been filled with the right unrolled index
print( look)
print(A[look])
}

\keyword{spatial}

