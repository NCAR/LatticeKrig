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

\name{LKrig.sim}

\alias{LKrig.sim}
\alias{LKrig.sim.conditional}
\alias{simConditionalDraw}
\title{
Functions for simulating a multi resolution process following the
Lattice Krig covariance model.

}
\description{
The fields are Gaussian and can be either simulated unconditionally or
conditional on the field values and a set of irregular locations.  }

\usage{
#
LKrig.sim(x1, LKinfo, M=1,just.coefficients = FALSE)
LKrig.sim.conditional( LKrigObj, M=1, x.grid= NULL, grid.list=NULL,
                           nx=80, ny=80,...,Z.grid=NULL, seed=42, verbose=FALSE)

simConditionalDraw(index = 1, LKrigObj, ghat, x.grid, Z.grid, PHIGrid,
                 seeds = 123, verbose = FALSE)

}
\arguments{


\item{grid.list}{Specifies a grid of spatial locations using the
\code{grid.list} format ( help(grid.list)). These are the locations
used to evaluate the fields generated from conditional simulation. The
default is to generate an 80X80 grid based on range of the
observations.}  \item{just.coefficients}{If TRUE just simulates the
coefficients from the Markov Random field.}


\item{LKinfo}{A list with components that give the information
describing a multi-resolution basis with a Markov random field used for
the covariance of the basis coefficients. This list is created in
\code{LKrig} and is returned as part of the output object or in a more
hands on manner directly using \code{LKrigSetup} (See section on
returned value below for this list's description.)}


\item{M}{Number of independent simulated fields.}

\item{nx}{Number of grid points in x coordinate for output grid.}

\item{ny}{Number of grid points in y coordinate for output grid.}

\item{LKrigObj}{An \code{LKrig} object, i.e. the output list returned by
\code{LKrig}.}

\item{seed}{Seed to set random number generator.}

\item{x1}{ A two column matrix of 2-dimension locations to evaluate
basis functions or the first set of locations to evaluate the
covariance function or the locations for the simulated process.  Rows
index the different locations: to be precise \code{x1[i,1:2]} are the
"x" and "y" coordinates for the \code{i} th location. }

\item{x.grid}{Locations to evaluate conditional fields. This is in the
form of a two column matrix where each row is a spatial location.}
	
\item{Z.grid}{The covariates that are associated with the x.grid
values. This is useful for conditional simulation where the fields are
evaluated at x.grid locations and using covariate values
Z.grid. Z.grid is matrix with columns indexing the different
covariates and rows indexed by the x.grid locations.}

\item{\dots}{Arguments to be passed to the LKrig function to specify
the spatial estimate. These are components in addition to what is in
the LKinfo list returned by LKrig.}

\item{verbose}{If TRUE prints out debugging information.}

\item{ghat}{The predicted surface at the grid.}

\item{index}{The index for the random seed to use in the vector \code{seeds}.}

\item{PHIGrid}{Basis function matrix at grid points.}

\item{seeds}{A vector of random seeds.}

} 

\details{
The simulation of the unconditional random field is done by
generating a draw from the multi-resolution coefficients using a
Cholesky decomposition and then multiplying by the basis functions to
evaluate the field at arbitrary points. Currently, there is no
provision to exploit the case when one wants to simulate the field on
a regular grid. The conditional distribution is a draw from the
multivariate normal for the random fields conditioned on the
observations and also conditioned on covariance model and covariance
parameters. If the nugget/measurement error variance is zero then any
draw from the conditional distribution will be equal to the
observations at the observation locations. In the past conditional
simulation was known to be notoriously compute intensive, but the major
numerical problems are finessed here by exploiting sparsity of the
coefficient precision matrix.
 
The conditional field is found using a simple trick based on the
linear statistics for the multivariate normal. One generates an
unconditional field that includes the field values at the
observations.  From this realization one forms a synthetic data set
and uses LKrig to predict the remaining field based on the synthetic
observations. The difference between the predicted field and the
realization (i.e. the true field) is a draw from the conditional
distribution with the right covariance matrix. Adding the conditional
mean to this result one obtains a draw from the full conditional
distribution.  This algorithm can also be interpreted as a variant on
the bootstrap to determine the estimator uncertainty. The fixed part
of the model is also handled correctly in this algorithm.  See the
commented source for \code{LKrig.sim.conditional} for the details of
this algorithm.

\strong{simConditionalDraw} is low level function that is called to generate each
ensemble member i.e. each draw from the conditional distribution. The large number of 
arguments is to avoid recomputing many common elements during the loop in generating these
draws. In particular passing the  basis function matrices avoid having to recompute the 
normalization at each step, often an intensive computation for a large grid. 


}

\value{

\strong{LKrig.sim:} A matrix with dimensions of \code{nrow(x1)} by
  \code{M} of simulated values at the locations \code{x1}.

\strong{LKrig.sim.conditional:} A list with the components. \describe{
\item{xgrid}{The locations where the simulated field(s) are
evaluated.}

\item{ghat}{The conditional mean at the xgrid locations.}

\item{g.draw}{A matrix with dimensions of \code{nrow(x.grid)} by
\code{M} with each column being an independent draw from the
conditional distribution.}




}
}
\author{
Doug Nychka
}

\seealso{
LKrig, mKrig, Krig, fastTps, Wendland
}
\examples{
# Load ozone data set
  data(ozone2)  
  x<-ozone2$lon.lat
  y<- ozone2$y[16,]
# Find location that are not 'NA'.
# (LKrig is not set up to handle missing observations.)
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
  LKinfo<- LKrigSetup( x,NC=20,nlevel=1, alpha=1, lambda= .3 , a.wght=5)
# BTW lambda is close to MLE 
# Simulating this  LKrig process
# simulate 4 realizations of process and plot them
# (these have unit marginal variance)
  xg<- make.surface.grid(list( x=seq( -87,-83,,40), y=seq(36.5, 44.5,,40)))
  out<- LKrig.sim(xg, LKinfo,M=4)
\dontrun{
  set.panel(2,2)
  for( k in 1:4){
    image.plot( as.surface( xg, out[,k]), axes=FALSE) }
}
  obj<- LKrig(x,y,LKinfo=LKinfo)
  O3.cond.sim<- LKrig.sim.conditional( obj, M=3,nx=40,ny=40) 
\dontrun{
  set.panel( 2,2)
  zr<- range( c(  O3.cond.sim$draw,  O3.cond.sim$ghat), na.rm=TRUE)
  coltab<- tim.colors()
  image.plot( as.surface( O3.cond.sim$x.grid, O3.cond.sim$ghat), zlim=zr)
  title("Conditional mean")
  US( add=TRUE)
  for( k in 1:3){
    image( as.surface( O3.cond.sim$x.grid, O3.cond.sim$g.draw[,k]),
              zlim=zr, col=coltab)
    points( obj$x, cex=.5)
    US( add=TRUE)
  }
  set.panel()
}

}

\keyword{spatial}


