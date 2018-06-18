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

\name{LKrig Internal}
\alias{LKrig.coef}
\alias{LKrig.lnPlike}
\alias{LKrig.traceA}
\alias{dfind2d}
\alias{dfind3d}
\title{
Internal functions for LatticeKrig package. 
}
\description{
Some internal functions for \code{LKrig} that estimate  the 
coefficients of the basis functions and compute the likelihood.
}
\usage{
LKrig.coef(Mc, wPHI, wT.matrix, wy, lambda, weights)
LKrig.lnPlike(Mc, Q, y, lambda, residuals, weights,
 sigma = NA, rho = NA)
LKrig.traceA( Mc, wPHI, wT.matrix, lambda, weights,NtrA, iseed=NA)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iseed}{Random seed used to generate the Monte Carlo
   samples. Keep the same to compare results with mKrig and also for
   multiple values of lambda.}
 
 \item{lambda}{The ratio of the nugget variance (sigma squared) to the
  parameter controlling the marginal variance of the process (called
  rho in fields).}

  \item{Mc}{Cholesky decomposition of regression matrix.}

  \item{NtrA}{Number of Monte Carlo samples to estimate trace. Default
  is 20 in LKrig.}

  \item{Q}{Precision matrix for coefficients.}

   \item{residuals}{Residuals from fitting spatial process.}

  \item{rho}{If lambda is not specified then the value of the marginal
  variance of the process.}

  \item{sigma}{If lambda is not specified then the values of the
  measurement error standard deviation.}

  \item{wPHI}{Weighted matrix of basis functions. See \code{LKrig}
  source for construction.}
  
  \item{wT.matrix}{Weighted matrix of fixed part of estimate.}

  \item{weights}{A vector that is proportional to the reciprocal
  variances of the errors. I.e. errors are assumed to be uncorrelated
  with variances sigma^2/weights. }

  \item{wy}{Weighted observations.}

  \item{y}{Spatial observation.}
}
\details{
\code{LKrig.coef} and \code{LKrig.lnPlike} are two low level functions
to find the basis function coefficients and to evaluate the
likelihood. The coefficents (\code{c.mKrig}) are also found because
they provide for shortcut formulas for the standard errors and MLE
estimates.  These coefficients are identical to the basis coefficients
(\code{c.coef}) found for usual Kriging in the mKrig
function. \code{LKrig.lnPlike} also finds the profile MLE of sigma and
rho given a fixed value for lambda (and \code{alpha} and
\code{a.wght}).  See the source for LKrig and also MLE.LKrig to see
how these functions are used.

\code{LKrig.traceA} finds an estimate of the effective degrees of
freedom of the smoothing matrix based a simple Monte Carlo scheme. The
smoothing matrix A is the matrix for fixed covariance parameters so
that y.hat = A y, where y.hat are the predicted values at the data
locations.  trace(A) is the effective degrees of freedom.  If e are
iid N(0,1) then the expected value of t(e)\% * \% A \% * \% e is equal
to the trace of A. This is the basis for estimating the trace and the
standard error for this estimate is based on \code{NtrA} independent
samples.
 
\code{dfind2d} is a fast FORTRAN subroutine to find nearest neighbors
within a fixed distance and is called by \code{Wendland.basis}. The
function \code{dfind3d} is currently not used but is intended for
future use to determine chordal distance between points on a sphere or
cylinder.

}
\value{
\code{LKrig.coef} a list with components d.coef the coefficients of
the spatial dirft and for covariates (Z) and c.coef the basis function
coefficients. The logical vector ind.drift from the LKrig object
indicates with components of d.coef are associated with the polynomial
spatial drift and which are other fixed spatial covariates.

\code{LKrig.lnPlike} has the components:
\describe{

\item{lnProfileLike:}{the log likelihood profiled for lambda, alpha
and a.wght}

\item{rho.MLE:}{the MLE of rho given lambda, alpha and a.wght}

\item{shat.MLE:}{ the MLE of sigma given lambda, alpha and a.wght}

\item{quad.form:}{the quadratic form in the exponent of the
 multivariate normal likelihood}

\item{lnDetCov:}{the log determinant of the covariance matrix in the
 likelihood}
}
}
\author{
Doug Nychka
}
\seealso{
LKrig, LKrig.basis
}
\keyword{spatial}
