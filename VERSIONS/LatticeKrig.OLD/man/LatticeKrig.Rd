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
% # Modified  NL Jan -> DWN Jan 2014

\name{LatticeKrig}
\alias{LatticeKrig}
\title{User-friendly spatial prediction and inference using a compactly supported
multi-resolution basis and a lattice model for the basis coefficients.
}
\description{ This is a simple version of the \code{LKrig} function
that creates a a complete \code{LKrig} object that can then used for
spatial prediction and inference. It uses the defaults from
\code{LKrig} and also uses \code{LKrigFindLambda} to choose a lambda
through likelihood maximization. Despite the simple syntax,
LatticeKrig still takes advantage of the multi-resolution feature of
\code{LKrig} and any \code{LKrig} parameter can be passed through the
function call.}

\usage{
LatticeKrig(x, y, Z=NULL, nu=1, nlevel=4, a.wght=4.01, ...)



}

\arguments{
  \item{a.wght}{The "a" spatial autoregressive parameter for a first
 order Markov Random field. This controls the spatial dependence and
 must be greater than or equal to 4. Default is a large spatial
 similar to thin plate spline type model. See Details.}

  \item{nlevel}{Number of levels for the multi-resolution basis. Each
level increases the number of basis functions by roughly a factor of
4.}
 
  \item{nu}{Smoothness of the multi-resolution process. If passed,
  alpha vector is proportional to \code{exp(-(1:nlevel)*nu)} and
  normalized to sum to one. Default is similar to a Whittle covariance
  (Matern with smoothness of 1).}

  \item{x}{Spatial locations of observations. Or the \code{LKrig}
  object for printing. }
 
  \item{y}{Spatial observations. No missing values are allowed.}
 
  \item{Z}{Linear covariates to be included in fixed part of the model
  that are distinct from the default first order polynomial in
  \code{x} (i.e. the spatial drift).}

  \item{\dots}{Additional arguments to pass to LKrig. }

}
\details{

This function is a wrapper that simples the model choices and calling
sequence to the more general function LKrig. It is useful for users
not fully familiar with the Lattice Krig methodology or those that
wish to try a default approach to get a quick look at a spatial
analysis. The default values are set to give about 4 times as many
basis functions as observations and to use four levels of
multi-resolution. The spatial correlation range is nearly stationatry
and set large to mimic a thin-plate spline. The smoothness mimic the
Whittle covariance function ( smoothness =1 for the Matern).  (see
LKrig.cov.plot to get a plot of the implied covariance function.)

Of course this function is pitched with all the caveats that
statistical model assumptions should always be checked and applying
generic methods to a specific problems without checking the
appropriateness can lead to misleading results. Details on the full
computations can be found in the \code{LKrig} man page.

The \code{lambda} parameter in \code{LKrig} is essential to the
Lattice Krig computation and an inappropriate value will result in
over or under fitting and incorrect interpolated values. The function
\code{LKrigFindLambda} is used within \code{LatticeKrig} to estimate a
lambda value from the data using maximum likelihood.

The main call inside \code{LatticeKrig} is to \code{LKrig} and so a
\code{LKrig} object is returned. Thus all of the functionality that
comes with \code{LKrig} objects such as \code{predict},
\code{summary}, \code{predictSurface}, etc. remain the same as
described in \code{LKrig}.  Also see the components \code{residuals}
and \code{fitted.values} in the returned object for these parts of the
fitted spatial model. The component \code{LKinfo} has all the details
that specfiy the basis functions and covariance model.

}
\author{
Doug Nychka
}

\seealso{
LKrig, LKrig.setup, LKrigFindLambda, LKinfo, LKrig.sim.conditional
}
\examples{
# Load ozone data set
  data(ozone2)  
  x<-ozone2$lon.lat
  y<- ozone2$y[16,]
# Omit locations that are not 'NA'.
# (LKrig is not set up to handle missing observations.)
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]

# thin plate spline-like model with the lambda parameter estimated by
# maximum likelihood. Default choices are made for a.wght, nlevel, NC
# and alpha.
\dontrun{
  obj<- LatticeKrig( x, y)
# summary of fit and a plot of fitted surface
  print( obj)
  surface( obj )
  US(add=TRUE)
  points(x)
# prediction standard errors
  out.se<- predictSE( obj, xnew= x)
}

###########################################################################
# Including a covariate (linear fixed part in spatial model)
##########################################################################

  data(COmonthlyMet)
  y.CO<- CO.tmin.MAM.climate
  good<- !is.na( y.CO)
  y.CO<-y.CO[good]
  x.CO<- as.matrix(CO.loc[good,])
  Z.CO<- CO.elev[good]

\dontrun{
  obj  <- LatticeKrig(x.CO,y.CO,Z=Z.CO)
  obj2 <- LatticeKrig(x.CO,y.CO)

# compare with and without linear covariates
  set.panel(1,2)
  surface(obj)
  US(add=TRUE)
  title("With Elevation Covariate")

  surface(obj2)
  US(add=TRUE)
  title("Without Elevation Covariate")

}
}
\keyword{spatial}

