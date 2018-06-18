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
% Earlier modifications DH 2/06/2012


\name{Radial.basis}
\alias{Radial.basis}
\alias{LKrig.cyl}
\alias{WendlandFunction}
\title{
Two dimensional radial basis functions based on a Wendland function.
}
\description{
Two dimensional radial basis functions based on a Wendland function and
using sparse matrix format to reduce the storage.
}
\usage{
Radial.basis(x1, centers, delta, max.points = NULL, mean.neighbor =
                 50, just.distance = FALSE, RadialBasisFunction =
                 "WendlandFunction",distance.type="Euclidean")
LKrig.cyl( U, R=1 ) 
WendlandFunction(d)
}

\arguments{

\item{x1}{A 2 column matrix of 2-d locations to evaluate the basis
functions. Each row of \code{x1} is a location.}

\item{centers}{A two column matrix specifying the basis function
centers.}

\item{d}{A vector of distances.}

\item{delta}{A vector of scale parameters for the basis functions. }

\item{max.points}{Maximum number of nonzero entries expected for the
returned matrix.}

\item{distance.type}{The distance metric. See the entry in
\code{LKrig} for details.}

\item{mean.neighbor}{Average number of centers that are within delta
of each x1 location. For centers on a regular grid this is often easy
to estimate.}

\item{just.distance}{For debugging, don't evaluate the RBF just return
the pairwise distance matrix.}

\item{R}{Radius used for cylinder coordinates.}

\item{RadialBasisFunction}{The positive definite kernel used to
generate the basis. The assumption is that this will take a
nonnegative argument and be zero outside [0,1].}

\item{U}{A two column matrix giving the locations in cylindrical
coordinates.  First column is the angle in degrees and second column is
the height.}

}
\details{

This function finds the pairwise distances between the points x1 and
centers and evaluates the function RadialBasisFunction at these
distances scaled by delta. In most applications delta is constant, but
a variable delta could be useful for lon/lat regular grids. The
Wendland function is for 2 dimensions and smoothness order 2. See
\code{WendlandFunction} for the polynomial form. This code has a very
similar function to the fields function \code{wendland.cov}.

In pseudo R code for delta a scalar \code{Wendland.basis} evaluates as

\preformatted{
  BigD<- rdist( x1,centers)
  WendlandFunction(BigD/delta)
}

The actual code uses a FORTRAN subroutine to search over distances
less than delta and also returns the matrix in sparse format.

The function \code{LKrig.cyl} transforms coordinates on a cylinder,
e.g. lon/lat when taken as a Mercator projection, and returns the 3-d
coordinates. It is these 3-d coordinates that are used to find distances
to define the radial basis functions. For points that are close this
"chordal" type distance will be close to the geodesic distance on a
cylinder but not identical.

}
\value{

For \code{Wendland.basis} a matrix in sparse format with number of
rows equal to nrow(x1) and columns equal to nrow(center).

}
\author{
Doug Nychka
}
\seealso{
LKrig.basis
}
\examples{
x<- cbind( runif(100), runif(100))
center<- expand.grid( seq( 0,1,,5), seq(0,1,,5))
# coerce to matrix
center<- as.matrix(center)

PHI<- Radial.basis(x, center, delta=.5)

# LKrig with a different radial basis function. 
# 
  data(ozone2)  
  x<-ozone2$lon.lat
  y<- ozone2$y[16,]
# Find location that are not 'NA'.
# (LKrig is not set up to handle missing observations.)
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
  obj<- LKrig(x,y,NC=30,nlevel=1, alpha=1, lambda=.01, a.wght=5)
  triweight<- function( d){
       ifelse( abs(d)<=1, (1-d^2)^3, 0)}   
  obj1<- LKrig(x,y,NC=30,nlevel=1, alpha=1, 
    lambda=.01, a.wght=5, RadialBasisFunction="triweight", overlap=1.8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

