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

\name{LKDist}
\alias{LKDist}
\alias{LKDistComponents}
\alias{LKDistGrid}
\alias{LKDistGridComponents}
\alias{LKGridCheck}
\alias{LKGridFindNmax}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find all  pairwise distances  within a 
maximum distance.
}
\description{
	These are the lower level functions to 
	compute the distances among two sets of locations but being limited to 
	distances less than a maximum threshold (see delta below ). These functions are useful for 
	generating a sparse matrix on distances and evaluating a compactly supported function (such as the Wendland). The location - location method supports the distance metrics:
	Euclidean, spherical, component-wise and Manhattan.
	 \code{LKDistComponent} and  \code{LKDistComponentGrid} return the coordinate-wise distances and are
	  useful for evaluating a tensor product basis functions. 
}
\usage{
LKDist(x1, x2, delta, max.points = NULL, mean.neighbor = 50,
             distance.type = "Euclidean")
LKDistComponents(x1, x2, delta, max.points = NULL, mean.neighbor = 50, 
    distance.type = "Euclidean")
    
LKDistGrid(x1, gridList, delta, max.points = NULL, mean.neighbor = NULL,
distance.type = "Euclidean", periodic)

LKDistGridComponents(
x1, gridList, delta,
  max.points = NULL, mean.neighbor = NULL,
  distance.type = "Euclidean") 

LKGridCheck(distance.type, x1, gridList )

LKGridFindNmax(n1, max.points, mean.neighbor, delta, gridList)
}

\arguments{
  \item{gridList}{A list with each component vector that specifies the
  	grid points for an equally spaced grid. Can have class
  	gridList. (See also help on gridlist).}	
  \item{n1}{Number of rows of x1.}
 % \item{x}{A two column matrix giving longitude/ latitude locations on a
 %	 sphere.}
  \item{x1}{
  	A matrix with rows indexing locations and columns indexing coordinates.

}
  \item{x2}{
  	A matrix with rows indexing locations and columns indexing coordinates.
}
  \item{delta}{
  	The maximum distance to find pairwise distances. 

}
  \item{max.points}{
  	Used for dynamically assigning matrix size this should be larger than  the total number of pairwise distances less than delta. 
}
  \item{mean.neighbor}{
  	Used for dynamically assigning matrix size this is the average number of 
  	points that are less that delta in distance to the \code{x1} locations.
}
\item{periodic}{A logical vector with length ncol( x1). If a 
	component is TRUE then that dimension is treated as periodic.}
  \item{distance.type}{
  	A text string either "Euclidean", "GreatCircle", "Chordal", "Manhattan".
}  	
%  \item{components}{ If TRUE the components of distance are returned for each %coordinate separately. (But all components must be less than \code{delta} )
%}
  
}
\value{
	\strong{LKDist} and \strong{LKDistGrid} a list representing a sparse matrix in spind format.
	
	\strong{LKDistComponent} and\strong{LKDistGridComponent} a list representing a sparse matrix  using the  spind format except the \code{ra} component is now a matrix.  The columns of ra 
	being the individual distances for each coordinate. 
	 
	\strong{directionCosine} a matrix where rows index the different points and 
	columns index x,y,z.
	
	\strong{LKGridFindNmax} returns the maximum number of nonzero elements expected in a pairwise distance matrix. 
	
	\strong{LKGridFindNmax} checks that the calling arguments are 
	compatible with the  pairwise distance computation. 
	
}

\author{
Doug Nychka
}

\seealso{
	\code{\link{spind2spam}}, \code{\link{spind2full}}
}
\examples{
set.seed( 123)	
x<- matrix( runif(100*2), 100,2)

DMatrix<- LKDist( x,x, delta=.1)
# coerce to spam matrix format
DMatrix2<- spind2spam( DMatrix)

# a grid
gridL<- list( x1= seq(0,1,.2), x2= seq( 0,2,.2) , x3= seq( -1,1,.2))
class(gridL)<- "gridList"	
x1<- cbind( runif( 100), runif(100)*2, 2*(runif( 100) -.5) )
look<- LKDistGrid( x1, gridL, delta=.45)
# check against rdist.
# look2<- rdist( x1, make.surface.grid(gridL))
# look2[ look2 >= .45] <- 0
# max( abs(look- look2)[look>0] )

# test of periodic option
 gridL<- structure(
            list( x1= seq(0,1,.02),
                  x2= seq( 0,1,.02)),
            class="gridList")
 look1<- LKDistGrid( rbind(c(0,0)), gridL, delta=.35,
                     periodic=c(TRUE,FALSE))
 look2<- spind2full(look1)
 image.plot( as.surface( gridL, look2) )
 
 look1<- LKDistGrid( rbind(c(0,0)), gridL, delta=.35,
                      periodic=c(TRUE,TRUE))
 look2<- spind2full(look1)
 image.plot( as.surface( gridL, look2) )  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

