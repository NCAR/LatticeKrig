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

\name{LKinfo}
\alias{LKinfo}
\alias{LKrig.setup}
\alias{LKrig.make.centers}
\alias{LKinfoUpdate}
\alias{LKrig.make.Normalization}
\alias{LKrig.make.a.wght}
\alias{LKrig.make.alpha}
\alias{LKrig.make.grid.info}
\alias{LKrigMRFDecomposition}
\title{
Specifying the Lattice Krig covariance object(\code{LKinfo}) and
related utility functions.   
}

\description{ The \code{LKinfo} object is a list that contains all the
details to specify an LKrig covariance model.  Typically this object
is created by supplying a few parameters to the function
\code{LKrig.setup}.  }

\usage{ LKrig.setup(x=NULL, NC=NULL, NC.buffer=5, nlevel,
grid.info=NULL, lambda=NA, sigma=NA, rho=NA, alpha=NA, nu=NULL,
a.wght=NA, overlap=2.5, normalize=TRUE, normalize.level=NULL,
edge=FALSE, rho.object=NULL, RadialBasisFunction="WendlandFunction",
distance.type = "Euclidean", V =diag(c(1, 1)), verbose=FALSE)

LKrig.make.centers( grid.info, nlevel, NC.buffer, distance.type)

LKrig.make.a.wght( a.wght, nlevel, mx, my)

LKrig.make.alpha( alpha, nu, nlevel)

LKrig.make.Normalization(mx,my, a.wght)

LKrigMRFDecomposition( mx,my,a.wght)

LKinfoUpdate(LKinfo, ... )

} \arguments{
 
\item{alpha}{A vector of length \code{nlevel} with the relative
variances for the different mulitresolution levels. }

\item{a.wght}{A vector of length \code{nlevel} that are the weights
for the central point in the Markov random field
specification. To be precise at level k the center point has weight
a.wght[k] with the 4 nearest neighbors given weight -1. The a.wght
must be greater than 4 for the fields to be stationary.}

\item{distance.type}{See entry in \code{LKrig}.}

\item{edge}{If FALSE no adjustments are made to edges.}

\item{grid.info}{A list with components \code{xmin, xmax, ymin, ymax,
 delta} that specifies the range and spacing for the coarsest level
 grid. NOTE: If \code{NC} is not specified \code{grid.info} must be
 passed.}

\item{lambda}{The "noise to signal ratio" or also known as the
smoothing parameter it is the parameter \code{lambda =
sigma^2/rho}. If specified then sigma and rho typically are estimated
in LKrig by maximum likelihood. If lambda is not specified then it is
set as \code{lambda = sigma^2/ rho}. Note that to evaluate the spatial
process model, e.g. using the function \code{LKrig.cov}, a value of
lambda is not needed and this argument can default to NA.}

\item{LKinfo}{An LKinfo object as described in the Details or that
created by \code{LKrig.setup}.}

\item{mx}{Number of lattice points in horizontal coordinate.}

\item{my}{Number of lattice points in vertical coordinate.}

\item{NC}{Maximum number of grid points in either the x or y
 directions to determine the basis function centers. If the region is
 square then this will be \code{NC} grid points in each dimension
 within the spatial domain defined by the data locations (or
 grid.info) giving a total \code{NC**2} grid points/basis functions
 for the coarsest level. There may be a larger number of total basis
 function based on how the buffer grid domain is set. }
 \item{NC.buffer}{Number of grid points to add beyond the spatial
 domain. This number is fixed at each level. } \item{nlevel}{Number of
 levels in multiresolution. Note that each subsequent level increases
 the number of basis functions within the spatial domain size by a
 factor of roughly 4. }

\item{overlap}{The amount of overlap between basis functions. If the
grid spacing in x and y is \code{delta} then the support of the basis
functions will be \code{overlap*delta}.  This value is hardwired to
2.5 in \code{LKrig} and also the default in \code{LKrig.setup}}
 
\item{normalize}{If TRUE the basis functions will be normalized to
give a marginal variance of one.}

\item{normalize.level}{A logical vector of length number of
multiresolution levels. If \code{normalize.level[i]==TRUE} then the
process at level i has its basis functions "normalized" to give a
constant marginal variance.  Default is that all values are TRUE. This
is an experimental argument to save on computation, which might be
used to avoid normalizing very fine levels that have minimal
artifacts, but many basis.}

\item{nu}{If passed alpha vector is defined as \code{
 exp(-2*(1:nlevel)*nu)} but then normalized to sum to one. }

\item{RadialBasisFunction}{ A character string with the name of the R
function used to generate the multiresolution basis. Default is
\code{WendlandFunction}, the Wendland covariance (order 2).}

\item{rho.object}{A prediction object to specify part of the marginal
variance for the process.  If omitted assumed to be the constant
one. Calling \code{predict(rho.object,x1)} should return a vector with
the values at the locations in x1. }


\item{rho}{ A scalar, the sill or marginal variance of the process.}

\item{sigma}{The measurement error standard deviation.}

\item{x}{Spatial locations for fitting surface. This is used to
determine ranges of the grid for basis functions so only two points
are required that bound the rest of the data locations. E.g.  x=
cbind( c( 0,1), c(0,1)) will set the domain to be the unit square.  }

\item{V}{See entry in \code{LKrig}.}

\item{verbose}{If TRUE print out intermediate information.}

\item{\dots}{Any additional arguments that can be used in LKrig.setup
argument}.} 

\details{The basis functions are two-dimensional radial
basis functions based on the Wendland covariance function, centered on
regular grid points and with the scaling tied to the spacing of the
grid points. These grids are created by the utility
\code{LKrig.make.centers} listed above.

For a basis at the coarsest level, the grid centers for the first
level are generated by expanding the sequences in the x and y
coordinates into a regular grid of center points and using the spacing
either determined from the domain ranges and \code{NC} or just from
the \code{delta} component in \code{grid.info}.  The following R code
demonstrates how to create a grid without any buffer grid points
(i.e. \code{NC.buffer} equal to zero) outside the spatial domain.

\preformatted{
   xgrid<- seq(grid.info$xmin,grid.info$xmax,grid.info$delta)
   ygrid<- seq(grid.info$ymin,grid.info$ymax,grid.info$delta)
   The.Lattice<- make.surface.grid( list( x=xgrid, y=ygrid))
} 

Note that the same spacing \code{delta} is used in both directions.

The subsequent levels use the same expressions but delta is reduced by
a factor 1/2 for each level.  The LKinfo object has several components
that help with the bookkeeping for the number of centers at each
level. 

The additional functions listed are typically called from the top
level function \code{LKrig.setup} and are used to complete the LKinfo
object.  \code{LKrig.make.centers} and \code{LKrig.make.grid.info}
setup the lattice and also the centers of the basis functions.

\code{LKrig.make.a.wght} and \code{LKrig.make.alpha} create the
covariance parameter lists based on defaults and what is passed.

\code{LKrig.make.Normalization} computes some matrices that speedup
normalizing the absis functions to have a constant marginal variance
function. These matrices depend on the a.wght parameter and so need to
be recomputed if a.wght is changed. The function \code{LKinfoUpdate}
is a more rigorous way to change just a few parameter in the LKinfo
object and has the advantage that when a.wght is changes the
normalization matrices are recalculated.

}

\value{

\strong{\code{LKrig.setup:}} Returns an \code{LKinfo} object.  This
 object is a "list of lists" and collects all the information to
 describe the covariance model and the observational model. In most case
 one would supply some of the important components of this and the
 functions \code{LKrig.setup} and \code{LKrig} will fill in the rest
 with their defaults. The basic design is to pass the LKinfo object to
 subsequent functions to avoid keeping track of many parts of the
 covariance specification. See for example \code{LKrig.cov}.

{Components of the \code{LKinfo} object}

\describe{

\item{mx my}{ Two vectors giving the size of the grids for each
level. At level \code{l} there are \code{ mx[l]*my[l]} basis
functions.}

\item{nlevel}{Number of multiresolution levels.}

\item{delta}{Spacing between grid point centers.}

\item{m}{Total number of basis functions.}

\item{offset}{Indices for coefficients by level. offset[k]+1 is the
position of the first coefficient in level k.}

\item{grid.info}{A list with min and max values for x and y grids and
the spacing.}

\item{grid}{ A list with nlevel components: grid[[k]] is a list with
components \code{x} and \code{y} specifying the grid for the basis
function centers at level k.}

\item{overlap}{Amount the basis functions overlap.}

\item{alpha}{A list with nlevel components with alpha[[k]] being the
scalar or matrix giving the values for alpha at the kth level.}

\item{a.wght}{A list with nlevel components with a.wght[[k]] being the
scalar or matrix giving the values for a.wght at the kth level.}

\item{lambda}{Value of lambda. A value of NA is OK if it has
not been estimated or set.}

\item{ sigma rho}{Values of the covariance parameters for measurement
standard deviation and the process marginal variance. These can be NA
if not yet determined.}

\item{normalize}{A logical value with TRUE being default to normalize the
basis function to have constant marginal variance.}

\item{ normalize.level}{ A logical vector of length nlevel to
determine which levels of basis functions should be normalized.}

\item{edge}{A logical parameter; if TRUE edge corrections are made to precision
matrix.}

\item{scalar.alpha}{Logical value with TRUE indicating that alpha is a
scalar at each level -- i.e. not a matrix.}

\item{scale.basis}{Logical if TRUE indicates rho is a prediction obj
defining a surface.}

\item{rho.object}{The prediction object that defines the value of rho
at every point in the domain according to the predict function.}

\item{RadialBasisFunction}{Name of radial basis function used as
template for generating basis functions.}  }


\strong{\code{LKrig.make.centers:}} Returns a list whose components
are a subset of those listed above:


}
\author{
Doug Nychka
}

\seealso{
LKrig, mKrig, Krig, fastTps, Wendland, LKrig.basis
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
  LKinfo<- LKrig.setup( x,NC=10,nlevel=2, alpha=c(1,.5),
                            lambda= NA , a.wght=c(5,5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

