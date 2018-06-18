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

\name{LKrig.fixed.component}
\alias{LKrig.fixed.component}

\title{
Creates fixed part of spatial model.
}
\description{Creates matrix of low order polynomial in the spatial
coordinates and adds any other spatial covariates that are part of
the linear model.}

\usage{ LKrig.fixed.component(x, Z = NULL, m,
distance.type="Euclidean") }

\arguments{

\item{x}{ A 2 column matrix of 2-d locations to evaluate the
   polynomial.}

\item{Z}{ A matrix specifying additional spatial covariates.} 

\item{m}{The order of the polynomial. Following the convention for
splines the polynomial will have maximum order (m-1). Throughout LKrig
m==2 is the default giving a linear polynomial.}

\item{distance.type}{The distance metric. See the entry in
\code{LKrig} for details.}

 
}
\details{
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
T.matrix<- LKrig.fixed.component(x, m=2)
# quadratic polynomial 
T.matrix<- LKrig.fixed.component(x, m=3)
}

\keyword{spatial}

