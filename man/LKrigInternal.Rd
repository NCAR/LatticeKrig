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

\name{LKrig Internal}
\alias{LKrig.coef}
\alias{LKrig.lnPlike}
\alias{LKrig.lnPlikeOLD}
\alias{LKrig.traceA}
\alias{dfind2d}
\alias{dfind3d}
\alias{LKrigUnrollZGrid}
\alias{LKrigMakewU}
\alias{LKrigMakewX}
\alias{createLKrigObject}

\title{
Internal functions for LatticeKrig package. 
}
\description{
Some internal functions for \code{\link{LKrig}} that estimate  the 
coefficients of the basis functions and compute the likelihood.
}
\usage{
createLKrigObject(x, y, weights = NULL, Z, X, U, LKinfo, verbose = FALSE)
LKrigMakewU(object, verbose = FALSE)
LKrigMakewX(object, verbose = FALSE)
LKrig.coef( GCholesky, wX, wU, wy, lambda,
collapseFixedEffect = FALSE,  verbose = FALSE)
LKrig.lnPlike( GCholesky, Q, quad.form, nObs, nReps, weights, LKinfo)
LKrig.lnPlikeOLD(Mc, Q, wy, residuals, weights, LKinfo)
LKrig.traceA(GCholesky, wX, wU, lambda, weights, NtrA, iseed = NA)
LKrigUnrollZGrid( grid.list, ZGrid=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{collapseFixedEffect}{If FALSE estimate fixed effects
            separately for each replicated data set. }
  \item{grid.list}{The grid for evaluating surface}
  \item{GCholesky}{SPAM cholesky decomposition of the "G" matrix.}

  \item{iseed}{Random seed used to generate the Monte Carlo
   samples. Keep the same to compare results with mKrig and also for
   multiple values of lambda.}
 
 \item{lambda}{The ratio of the nugget variance (sigma squared) to the
  parameter controlling the marginal variance of the process (called
  rho in fields).}
 \item{LKinfo}{The LKinfo object. See help(LKinfo)}
%  \item{m}{Parameter specifying the degree of polynomial is (m-1).}

  \item{Mc}{Cholesky decomposition of regression matrix.}

  \item{NtrA}{Number of Monte Carlo samples to estimate trace. Default
  is 20 in LKrig.}
  
  \item{nObs}{Number of observations.}
  
  \item{nReps}{Number of replicate fields. }
  
  \item{object}{The LKrig object.}

  \item{Q}{Precision matrix for coefficients.}
  \item{quad.form}{The part of the log likelihood that is a quadratic form.
  (This is typically found in \code{LKrig.coef}.)}

   \item{residuals}{Residuals from fitting spatial process.}
   
\item{U}{The matrix that maps the d.coef coefficients of the fixed component
  	(typically a low order polynomial) part of the observation model.}
  	
 	
\item{verbose}{If TRUE intermediate debugging information is printed.}
   
\item{weights}{A vector that is proportional to the reciprocal
  variances of the errors. I.e. errors are assumed to be uncorrelated
  with variances sigma^2/weights. }
  
\item{wU}{Weighted U matrix the fixed part of the model.}
\item{wX}{Weighted X matrix (in spam format) related to nonparametric (stochastic) part of 
	model. Here weights refer to the sqrt(weights).
	
	NOTE: predicted values are  U\%*\%d.coef + X\%*\%c.coef
	
	}
	
  \item{wy}{Weighted observations.}
  
  
  \item{X}{The matrix that maps the c.coef coefficients into the nonparametric component
  	(spatial process) part of the observation model.}

 \item{x}{Matrix of spatial locations passed to LKrig.}
 \item{y}{Vector or matrix of observations passed to LKrig.}
 \item{Z}{A matrix of covariates.}

  \item{ZGrid}{A list or array with the covariates on the same grid as that specified by 
the  grid.list argument. }

}
\details{

The LatticeKrig article can be used as a reference for the matrix computations
and the G matrix from those formulas figures prominently. The GCholesky object
in these functions is the cholesky decompoistion of this matrix. For
compatibility with older version of this package this object may also be named
as \code{Mc} ( Cholesky of the M  matrix) but the user should not identify this M with that in the article. Ideally all coding using Mc should be changed to GCholesky.

\code{createLKrigObject} Based on the arguments passed into LKrig forms the 
prototype LKrig object. This object is added to as one computes additional steps in the LKrig function. 
	
\code{LKrigMakewU} and \code{LKrigMakewX}  construct the weighted U and X matrices from what is passed. In the case of observations that are point locations wU is found the weights and using the fixedFunction and wX is found from the weights and the multiresolution basis functions. Note that X and wX are assumed to be in spam 
sparse matrix format. 

\code{LKrig.coef} and \code{LKrig.lnPlike} are two low level functions
to find the basis function coefficients and to evaluate the
likelihood. The coefficients (\code{c.mKrig}) are also found because
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

\code{LKrigDefaultFixedFunction} Is called to construct the fixed part of the 
spatial model. The default is a polynomial of degree (m-1).
}
\value{
\describe{
\item{LKrig.coef}{a list with components d.coef the coefficients of
the spatial dirft and for covariates (Z) and c.coef the basis function
coefficients. The logical vector ind.drift from the LKrig object
indicates with components of d.coef are associated with the polynomial
spatial drift and which are other fixed spatial covariates.}

\item{LKrig.lnPlike}{has the components:}


\item{lnProfileLike}{the log likelihood profiled for lambda, alpha
and a.wght}

\item{rho.MLE}{the MLE of rho given lambda, alpha and a.wght}

\item{shat.MLE}{ the MLE of sigma given lambda, alpha and a.wght}

\item{quad.form}{the quadratic form in the exponent of the
 multivariate normal likelihood}

\item{lnDetCov}{the log determinant of the covariance matrix in the
 likelihood}


\item{LKrigDefaultFixedFunction}{A matrix with dimension nrow(x) and columns
of the number of polynomial terms  and the number of columns of Z if given.}
}
}
\references{
Nychka, D., Bandyopadhyay, S., Hammerling, D., Lindgren, F., & Sain, S. (2015). A multiresolution gaussian process model for the analysis of large spatial datasets.{Journal of Computational and Graphical Statistics}, 24(2), 579-599.
}
\author{
Doug Nychka
}
\seealso{
LKrig, LKrig.basis
}
\keyword{spatial}

