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
% Earlier modifications DH 2/06/2012

\name{LKrigDefaultFixedFunction}
\alias{LKrigDefaultFixedFunction}
\alias{LKrigPeriodicFixedFunction}
\alias{predictLKrigFixedFunction}
\title{
Creates fixed part of spatial model.
}
\description{Creates matrix of low order polynomial in the spatial
coordinates and adds any other spatial covariates that are part of
the linear model.}

\usage{
 LKrigDefaultFixedFunction(x, Z = NULL, m=2,distance.type="Euclidean") 
 LKrigPeriodicFixedFunction(x, Z = NULL, m = 2, distance.type = "Euclidean")
 
predictLKrigFixedFunction(object, xnew, Znew = NULL, drop.Z = FALSE,
collapseFixedEffect = FALSE)
}

\arguments{

\item{collapseFixedEffect}{If FALSE the fixed part of the model is found seperately for each replicated data set. If TRUE the estimate is polled across replicates.}
 
\item{drop.Z}{If TRUE only spatial drift is evaluated the contribution for covariates
is omitted.}

\item{distance.type}{The distance metric. See the entry in
\code{LKrig} for details.}
 
\item{m}{The order of the polynomial. Following the convention for
splines the polynomial will have maximum order (m-1). Throughout LKrig
m==2 is the default giving a linear polynomial.}

\item{object}{A LKrig object.}

\item{x}{ A 2 column matrix of 2-d locations to evaluate the
   polynomial.}

\item{xnew}{Locations for predictions.}

\item{Z}{ A matrix specifying additional spatial covariates.} 

\item{Znew}{Same as Z.}

}
\details{

\strong{LKrigDefaultFixedFunction} This function creates the regression matrix for the 
fixed part of the spatial model. The default is a low order polynomial regression matrix 
of degree m-1. To this matrix are bound as columns any covariates passed as Z. 
Typically one would not need to modify this function. 
For more exotic fixed part models one can specify create and then specify a different 
function. See LKrig.setup and LKrig. NOTE: If the argument for this function is passed as NULL
then the subsequent computations do not include a fixed part in the model.

\strong{LKrigDefaultFixedFunction} This is same as LKrigDefaultFixedFunction except the first coordinate is ignored. i.e. it is assumed to be periodic so adding a polynomial does not make sense.

\strong{predictLKrigFixedFunction}
This function is simple, but is introduced to make the code modular and
to handle the case for cylindrical geometry where only latitude should
have a spatial term (to preserve periodicity in longitude).
}

\value{
A matrix where rows index the locations and columns are the different
spatial polynomial and covariates. 
}

\author{
Doug Nychka
}

\seealso{
LKrig.basis, LKrig
}
\examples{
x<- matrix( runif(100), nrow=50)
# linear polynomial 
T.matrix<- LKrigDefaultFixedFunction(x, m=2)
# quadratic polynomial 
T.matrix<- LKrigDefaultFixedFunction(x, m=3)
}

\keyword{spatial}


