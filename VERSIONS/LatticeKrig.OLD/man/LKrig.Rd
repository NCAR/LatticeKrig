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

\name{LKrig}
\alias{LKrig}
\alias{predict.LKrig}
\alias{predictSE.LKrig}
\alias{print.LKrig}
\alias{print.LKinfo}
\alias{surface.LKrig}
\alias{predictSurface.LKrig}
\title{Spatial prediction and inference using a compactly supported
multi-resolution basis and a lattice model for the basis coefficients.
}
\description{ A variation of Kriging with fixed basis functions that
uses a compactly supported covariance to create a regular set of basis
functions on a grid. The coefficients of these basis functions are
modeled as a Gaussian Markov random field (GMRF). Although the
precision matrix of the GMRF will be sparse the model can still
exhibit longer ranges of spatial dependence. The multi-resolution
feature of this model allows for the approximation of a wide variety
of standard covariance functions.  }

\usage{

LKrig(x, y = NULL, weights = rep(1, nrow(x)), Z = NULL,
                 LKinfo = NULL, iseed = 123, NtrA = 20, use.cholesky =
                 NULL, return.cholesky = TRUE, wPHI = NULL, return.wPHI
                 = TRUE, NC, nlevel, a.wght, alpha, nu = NULL, lambda =
                 NA, sigma = NA, rho = NA, rho.object = NULL, overlap =
                 2.5, normalize = TRUE, edge = FALSE,
                 RadialBasisFunction = "WendlandFunction", V =
                 diag(c(1, 1)), distance.type = "Euclidean", verbose =
                 FALSE)
\method{predict}{LKrig}( object, xnew=NULL,Znew = NULL, drop.Z=FALSE,
return.levels = FALSE, ...)

\method{predictSE}{LKrig}( object, xnew=NULL, Znew = NULL, drop.Z =
FALSE, verbose = FALSE, ...)

\method{surface}{LKrig}( object, ...)

\method{predictSurface}{LKrig}( object, ...)

\method{print}{LKrig}( x, digits=4, ...)

}

\arguments{
  
  \item{alpha}{ Variance weights for each level of resolution. This
     can be a vector or a list with \code{nlevel} components. As a
     list each component is a matrix specifying a weight for each
     coefficient.  See Details.}

  \item{a.wght}{The "a" spatial autoregressive parameter for a first
 order Markov Random field. This controls the spatial dependence and
 must be greater than or equal to 4. For \code{a.wght} = 4 normalize
 should be FALSE.  If there are multiple levels this can be a
 vector. See Details.}


  \item{distance.type}{The type of distance metric (as a text string)
     used in determining the radial basis functions.  The current
     choices are "Euclidean" and "cylinder" with cylinder being used
     to have a periodic boundary condition in the x (horizontal)
     coordinate.}  \item{digits}{Number of digits in printed output.}
   
  \item{drop.Z}{If true the fixed part will only be evaluated at the
spatial drift polynomial part of the fixed model. The contribution
from the other, Z, covariates in the fixed component will be omitted.}

  \item{edge}{If TRUE an adjustment is made in the GMRF weights for
  edge effects.}

 \item{iseed}{Random seed used in the Monte Carlo technique for
  approximating the effective degrees of freedom (trace of the
  smoothing matrix). If NA, no seed is set. }

 \item{lambda}{The ratio of the nugget variance (called sigma squared
  in \code{fields} and \code{LatticeKrig}) to the parameter
  controlling the marginal variance of the process (called rho in
  \code{fields} and \code{LatticeKrig}). If sigma and rho are both
  specified then \code{lambda == sigma^2/rho}.}

 \item{LKinfo}{An object whose components specify the LatticeKrig
  spatial model.  This is usually created by the function
  \code{LKrig.setup}. If NULL, this object is created and returned as
  a component of the LKrig object. }

   \item{NC}{Controls the maximum number of grid points. For the first
 level there will be NC grid points in the larger dimension of the
 rectangular spatial domain. If the domain is square there will be
 NC*NC basis functions/lattice points in the coarsest level of
 resolution. If domain is rectangular, the smaller dimension will have
 less than NC points.}

  \item{nlevel}{Number of levels for the multiresolution basis. Each
level increases the number of basis functions by roughly a factor of
4.}

  \item{normalize}{If \code{TRUE} basis functions are normalized so
    that the marginal variance of the process covariance is constant
    and equal to rho. This normalization avoids some of the edge and
    periodic artifacts from using a discrete set of basis functions.}

  \item{NtrA}{Number of random samples used in Monte Carlo method for
determining effective degrees of freedom. }
 
 \item{nu}{If passed, alpha vector is proportional to \code{
exp(-(1:nlevel)*nu)} and normalized to sum to one.}

\item{object}{The \code{LKrig} object.}

\item{overlap}{The overlap between basis functions. This scaling is
    based on centers being spaced 1 unit distance apart and the
    Wendland function decreasing to zero at 1 unit distance. A
    scaling/overlap of 2.5 (the default) implies that the support of
    the basis functions will extend to a disc of radius 2.5.  We
    recommend that this parameter not be changed unless one is sure of
    the effect on the implied spatial covariance.}

\item{return.cholesky}{If TRUE the Cholesky decomposition is included
    in the output list (with the name \code{Mc}). This is needed for
    some of the subsquent computations such as finding prediction
    standard errors.  Set this argument to FALSE to avoid much larger
    objects when the decomposition is not needed. This option is often
    paired with a subsequent call to \code{LKrig} with
    \code{use.cholesky}. See the MLE.LKrig source code for details. }

\item{return.levels}{If TRUE the predicted values for each level are
returned as columns in a matrix. }

  \item{RadialBasisFunction}{An R function that defines the radial
            basis function to generate the multiresolution basis. Is
            assumed that this function is compactly supported on the
            interval [0,1]. See WendlandFunction for an example.}

  \item{rho}{Value of rho (sill variance) to use for spatial
  estimate. If omitted this is estimated as the MLE based on the value
  of lambda.}

  \item{rho.object}{A prediction object that defines a spatially
           varying component for the marginal variance of the
           process. The object should be such that the code
           \code{predict(rho.object, x1)} will evaluate the function
           at 2-d locations \code{x1}. See Details below.  }

  \item{sigma}{Value of sigma (nugget variance) to use for the spatial
            estimate.  If omitted this is estimated as the MLE based
            on the value of lambda.}  \item{wPHI}{If not NULL then use
            as the "wPHI" regression matrix in the computations. This
            is only used when just lambda is being changed from a
            previous call.}

  \item{return.wPHI}{If TRUE return the wPHI regression matrix. This
  is most likely to be used for a subsequent call to LKrig that avoids
  recreating this matrix.}

  \item{verbose}{If \code{TRUE} print out intermediate results and
  messages.}
  

  \item{use.cholesky}{Use the symbolic part of the Cholesky
  decomposition passed in this argument.}

  \item{V}{A matrix used to scale and rotate coordinates to give an
anisotropic basis function and covariance model. V is the 2X2 matrix
describing the inverse linear transformation of the coordinates before
distances are found.  In R code this transformation is: \code{x1 \%*\%
t(solve(V))} Default is NULL, that is the transformation is just
dividing distance by the scalar value \code{theta}. See Details below.
If one has a vector of "theta's" that are the scaling for each
coordinate then just express this as \code{V = diag(theta)} in the
call to this function.}

  \item{weights}{A vector that is proportional to the reciprocal
  variances of the errors. I.e. errors are assumed to be uncorrelated
  with variances sigma^2/weights. }

  \item{xnew}{Matrix of locations for prediction.}
 
  \item{x}{Spatial locations of observations. Or the \code{LKrig}
  object for printing. }
 
  \item{y}{Spatial observations.}
 
  \item{Z}{Linear covariates to be included in fixed part of the model
  that are distinct from the default first order polynomial in
  \code{x} (i.e. the spatial drift).}

  \item{Znew}{Values of covariates, distinct from the spatial drift
  for predictions of data locations.}

  \item{\dots}{Additional arguments to pass to generic methods or from
 LatticeKrig to LKrig.  See help(predictSurface) for the details in
 calling LKrig.predictSurface}.
}
\details{
This method combines compactly supported basis
functions and a Markov random field covariance model to provide
spatial analysis for large data sets. The model for the spatial
field (or spatial process) is
 
 f(x) = N(x) + Z d + sum Phi.j(x) c.j.

x is a location in two dimensions, N(x) is a low order (linear)
polynomial, Z is a matrix of spatial covariates and d a coefficient
vector. Phi.j for 1 <= j <= m is a set of fixed basis functions and
c.j the coefficients. The variance of f(x) is given by the parameter
\code{rho} throughout this package. As explained below the process f
is a sum of \code{nlevel} independent processes that have different
scales of spatial dependence.  The \code{alpha} gives the relative
weighting between these processes.  Thus, the minimum set of
parameters needed to describe the covariance of f are the integer
\code{NC}, two scalars \code{rho} and \code{a.wght} and a vector
\code{alpha}. \code{alpha} has length the number of multiresolution
levels but we recommend that it be constrained sum to one.  Thus in
total there are a 1 + 2 + (nlevel - 1) parameters for a minimal
specification of the covariance.  Note that this parsimonious
specification results in a covariance that is close to being
stationary and isotropic when \code{normalize} is \code{TRUE}.  An
additional constraint on \code{alpha} is to make the weights
\code{alpha[j]} proportional to \code{exp( - 2*j*nu)} where \code{nu}
controls decay of the alpha weights. There is some theory to suggest
that \code{nu} is analgous to the smoothness parameter from the Matern
family (e.g. nu=.5 approximates the exponential). In this case the
covariance model requires just four parameters, \code{NC, rho, a.wght,
nu}.

The data model is

 Y.k = f(x.k) + e.k

Y.k are (scalar) observations made at spatial locations x.k with e.k
uncorrelated normal errors with variance sigma^2/weights.  Thus there
is a minimum of one new parameter from the data model: sigma. Note
that prediction only depends on the ratio \code{lambda = sigma^2/ rho}
and not surprisingly lambda plays a key role in specifying and fitting
a spatial model. Also taken with the model for f the minimum
parameters needed for a spatial prediction are still four \code{NC,
a.wght, nu} and \code{lambda}. For fixed lambda there are closed form 
expressions for the MLEs for sigma and rho. LKrig exploits this feature
by depending on lambda and then computing the MLEs for sigma and rho. 

\strong{Spatial prediction:} The basis functions are assumed to be
fixed and the coefficients follow a multivariate Gaussian
distribution. Given this spatial model for f, it is possible to
compute the conditional expectation of f given Y and also maximize the
likelihood for the model parameters, lambda, alpha, and a.wght.  This
setting is known as fixed rank Kriging and is a common strategy for
formulating a spatial model. Typically fixed rank Kriging is used to
reduce the dimension of the problem by limiting the number of basis
functions. We take a different approach in allowing for models that
might even have more basis functions than observations. This provides
a spatial model that can come close to interpolating the observations
and the spatial process is represented with many degrees of freedom.
The key is to make sure the model ingredients result in sparse
matrices where linear algebra is required for the computations. By
doing so in this package it is possible to compute the estimates and
likelihood for large numbers of spatial locations. This model has an
advantage over covariance tapering or compactly supported covariance
functions (e.g. \code{fastTps} from \code{fields}), because the
implied covariance functions can have longer range correlations.

\strong{Radial basis functions ( Phi.j ) :}
 The basis functions are two-dimensional radial basis functions (RBFs)
 that are derived from scaling and translations of a single covariance
 function. The default in LatticeKrig is to use the Wendland compactly
 supported stationary covariance (order 2 for 2 dimensions) that is
 scaled to be zero beyond a distance of 1.  For d the distance between
 spatial locations, this Wendland function has the standard form:

\code{(1 - d)^6 * (35 * d^2 + 18 * d + 3))/3} for d in [0,1]

  0  otherwise.

 For a single level the RBFs are
centered at a regular grid of points and with radial support
\code{delta*overlap} where \code{delta} is the spacing between grid
points. We will also refer to this grid of centers as a lattice and
the centers are also referred to as "nodes" in the RBF literature. The
overlap for the Wendland has the default value of 2.5 and this
represents a compromise between the number of nonzero matrix elements
for computation and the shape of the covariance functions.

To create a multi-resolution basis, 
each subsequent level is based on a grid with delta divided by
2. See the example below and \code{help(LKrig.basis)} for more details.
For multiple levels the basis functions can be grouped according to
the resolution levels and the coefficients can be grouped in a similar
manner. 
There is one important difference in the basis construction -- a
normalization -- and this aspect makes it different from a simple
radial basis function specification and is described below.

\strong{Markov random field (GMRF) for the coefficients (c.j) :}
Because the coefficients are identified with locations on a
lattice  it is easy to formulate a Markov random field for their
distribution based on the relationship to neighboring lattice points.
The distribution on the basis function coefficients is a multivariate
 normal, with a mean of zero and the the precision matrix, Q, (inverse
 of Q is the covariance matrix).  Q is partitioned in a block diagonal
 format corresponding to the organization of the basis functions and
 coefficients into levels of resolution. Moreover, coefficients at
 different levels are assumed to be independent and so Q will be block
 diagonal. If \code{nlevels} are specified, the ith block has a
 precision matrix based on a spatial autoregression with
 \code{a.wght[i]} being related to the spatial autoregressive parameter(s).
 Schematically in the simplest case the weighting for an interior
 lattice point and its four neighbors is
\preformatted{ 
.   .    .      .     .                                       

.   .   -1      .     . 
                     
.  -1  a.wght  -1     .
                      
.   .   -1      .     .  
  
.   .    .      .     .
Interior point      
}

The fundamental idea is that these weights applied to each point in
the lattice will result in a lattice of random variables that are
independent.  The specific precision matrix for each block (level),
Q.i, is implemented by \code{LKrig.MRF.precision}.  In the case when
alpha is a scalar, let C.i be the vector of basis coefficents at the
ith level then we assume that \code{ B \%*\% C.i} will be independent
N(0,1) random variables. By elementary properties of the covariance it
follows that the precision matrix for C.i is Q.i= t(B)\%*\%B. Thus,
given B one can determine the precision matrix and hence the
covariance matrix.  Each row of B, corresponding to a point in the
lattice in the interior, is "a" (\code{a.wght[i]}) on the diagonal and
-1 at each of the four nearest neighbors of the lattice points. Points
on the edges and corners just have less neighbors but get the same
weights.
 
This description is a spatial autoregressive model (SAR). The matrix Q
will of course have more nonzero values than B and the entries in Q
can be identified as the weights for a conditional autoregressive
model (CAR). Moreover, the CAR specification defines the neighborhood
such that the Markov property holds.  Values for \code{a.wght[i]} that
are greater than 4 give reasonable covariance models.  Moreover
setting a.wght[i] to 4 and \code{normalize} to FALSE in the call to
LKrig will give a thin-plate spline type model that does not have a
range parameter. An approximate strategy, however, is to set a.wght
close to 4 and get some benefit from the normalization to reduce edge
effects.

\strong{Multiresolution process} Given basis functions and
 coefficients at each level we have defined a spatial process g.i that
 can be evaluated at any location in the domain. These processes are
 weighted by the parameter vector alpha and then added together to
 give the full process.  It is also assumed that the coefficients at
 different resolution levels are independent and so the processes at
 each level are also independent. The block diagonal structure for Q
 does not appear to limit how well this model can approximate standard
 spatial models and simplifies the computations.  If the each g.i is
 normalized to have a marginal variance of one then g will have a
 variance that is the sum of the alpha parameters. Usually it is
 useful to constrain the alpha parameters to sum to one and then
 include an additional variance parameter, rho, be the marginal
 variance for g.  So the full model for the spatial process used in
 LatticeKrig is

g(x) = sqrt(rho) * sum.i sqrt(alpha[i]) * g.i(x)

The specification of the basis and GMRF is through the components of
the object LKinfo, a required component for many LatticeKrig
functions. High level functions such as \code{LKrig} only require a
minimal amount of information and combined with default choices create
the LKinfo list. One direct way to create the complete list is to use
\code{LKrig.setup} as in the example below.  For \code{nlevel==1} one
needs to specify \code{a.wght}, \code{NC}, and also \code{lambda}
related to the measurement error variance.  For a multiresolution
setup, besides these parameters, one should consider different values
of \code{alpha} keeping in mind that if this vector is not constrained
in some way ( e.g. \code{sum(alpha)==1}) it will not be identifiable
from \code{lambda}.


 \strong{The covariance derived from the GMRF and basis functions:} An
important part of this method is that the GMRF describes the
coefficients of the basis functions rather than the field itself. Thus
in order to determine the covariance for the observed data one needs
to take in account both the GMRF covariance and the expansion using
the basis functions. The reader is referred to the function
\code{LKrig.cov} for an explicit code that computes the implied
covariance function for the process f.  Formally, if P is the
covariance matrix (the inverse of Q) for the coefficients then the
covariance between the field at two locations x1 and x2, will be

   sum_ij P_ij Phi.i(x1) Phi.j(x2)

 Moreover under the assumption that coefficients at different levels
are independent this sum can be decomposed into sums over the separate
levels.  The function \code{LKrig.cov} evaluates this formula based on
the LKrig object (\code{LKinfo}) at arbitrary groups of locations
returning a cross covariance matrix. \code{LKrig.cov.plot} is a handy
function for evaluating the covariance in the x and y directions to
examine its shape.  The function \code{LKrig.cov} is also in the form
to be used with conventional Kriging codes in the fields package
(loaded by LatticeKrig) \code{mKrig} or \code{Krig} and can be used
for checking and examining the implied covariance function.

\strong{Normalizing the basis functions} The unnormalized basis
functions result in a covariance that has some non-stationary
artifacts (see example below). For a covariance matrix P and for any
location x one can evaluate the marginal variance of the process using
unnormalized basis functions for each multiresolution level.  Based
this computation there is a weighting function, say w.i(x), so that
when the normalized basis w.i(x) Phi.i(x) is used the marginal
variance for the multiresolution process at each level will be one.
This makes the basis functions dependent on the choice of Q and
results in some extra overhead for computation.  But we believe it is
useful to avoid obvious artifacts resulting from the use of a finite
spatial domain (edge effects) and the discretization used in the basis
function model. This is an explicit way to make the model stationary
in the marginal variance with the result that the covariance also
tends to be closer to a stationary model. In this way the
discretization and edges effects associated with the GMRF can be
significantly diminished.

The default in \code{LKrig} is \code{normalize = TRUE}. It is an open
question as to whether all levels of the multi-resolution need this
kind of explicit normalization. There is the opportunity within the
\code{LKrig.basis} function to only normalize specific levels with the
\code{normalize} being extended from a single logical to a vector of
logicals for the different levels.  To examine what these edge effect
artifacts are like the marginal variances for a 6 by 6 basis is
included at the end of the Examples Section.

\strong{Nonstationary and anisotropic modifications to the covariance}
Given that the process at each level has been normalized to have
marginal variance of one there are several other points where the
variance can be modified. The variance at level i is scaled by the
parameters \code{alpha[i]} and the marginal variance of the process is
scaled by \code{rho}.  Each of these can been extended to have some
spatial variation and thus provide a model for nonstationarity.

An option in specifying the marginal variance is to prescribe a
spatially varying multiplier. This component is specified by the
object \code{rho.object}. By default this is not included (or assumed
to be identically one) but, if used, the full specification for the
marginal variance of the spatial process at location \code{x} is
formally: \code{rho * predict(rho.object, x) * sum( alpha)} There is
then a problem of identifiability between these and a good choice is
to constrain \code{sum(alpha) ==1} so that \code{rho*
predict(rho.object, x)} is associated with the marginal variance of
the full spatial process.

A second option is to allow the alpha variance component parameters to
vary across the lattices at each level. For this case alpha is a list
with \code{nlevel} components and each component being a matrix with
the same dimensions as the lattice at that level.  The SAR weight
matrix is taken to be the usual weights but each row is scaled by the
corresponding value in the alpha weight matrix. To be precise the
weight matrix is given in psuedo R code by

\code{diag( 1/sqrt( c(alpha[[i]]))\%*\%B}

leading to the precision matrix at level i of

  Q.i = t(B) \%*\% diag( c(1/alpha[[i]]))\%*\%B

In this code note that the matrix of weights, \code{alpha[[i]]} is
being stacked as a larger vector to match the implicit indexing of the
basis coefficients.

LKrig also has the flexibility to handle more general weights in the
GMRF.  This is accomplished by \code{a.wght} being a list with as many
components as levels. If each component is a vector of length nine
then these are interpreted as the weights to be applied for the
lattice point and its 8 nearest neighbors (see help(
LKrig.MRF.precision) and the commented source code for details). If
each component is a matrix then these are interpreted as the
(nonstationary) center weights for each lattice point. Finally if the
component is an array with three dimensions this specifies the center
and 8 nearest neighhors for every point in the lattice.  At this point
the choice of these weights beyond a stationary model is experimental
and we will defer further documentation of these features to a future
version. 

\strong{The smoothing parameter lambda and effective degrees of
freedom} Consistent with other fields package functions, the two main
parameters in the model, sigma^2 and rho are parameterized as lambda =
sigma^2/rho and rho. The MLEs for rho and sigma can be written in
closed form as a function of lambda and these estimates can be
substituted into the full likelihood to give a concentrated version
that can numerically be maximized over lambda. The smoothing parameter
lambda is best varied on a log scale and is sometimes difficult to
interpret independent of the particular set of locations, sample size
and covariance. A more useful interpretation of lambda is through the
effective degrees of freedom and an approximate value is found by
default using a Monte Carlo technique. The effective degrees of
freedom will vary with the dimension of the fixed regression part of
the spatial model ( typical 3 = constant + linear terms) and the total
number of observations. It can be interpreted as the approximate
number of "parameters" needed to represent the spatial prediction
surface. For a fixed covariance model the spatial estimate at the
observation locations can be represented as f hat = A(lambda) y where
A is a matrix that does not depend on observations.  The effective
number of degrees of freedom is the trace of A(lambda) in analogy to
the least squares regression "hat" matrix and has an inverse,
monotonic relationship to lambda. The Monte Carlo method uses the fact
that if e are iid N(0,1) E( t(e) A(lambda) e) = trace( A(lambda)).

\strong{Descriptions of specific functions and objects:}

\strong{\code{LKrig}: } Find spatial process estimate for fixed
covariance specificed by \code{ nlevel, alpha, a.wght, NC, and lambda}
or this information in an \code{LKinfo} list.

\strong{\code{predict.LKrig, predictSE.LKrig}: } 
These functions evaluate the model at the the data locations or at
\code{xnew} if it is included. Note the use of the \code{drop.Z}
argument to either include the covariates or just restrict the
computation to the spatial drift and the smooth component. If
\code{drop.Z} is FALSE then then \code{Znew} needs to be included for
predictions off of the observation locations. The standard errors are
computed under the assumption that the covariance is known, that it is
the TRUE covariance for the process, and both the process and
measurement errors are multivariate normal. The formula that is used
involves some recondite shortcuts for efficiency but has been checked
against the standard errors found from an alternative formula in the
fields Krig function. (See the script Lkrig.se.tests.R in the tests
subdirectory for details.)


}
\value{

\strong{\code{LKrig: }}  
An LKrig class object with components for evaluating the estimate at
arbitrary locations, describing the fit and as an option (with
\code{Mc.return=TRUE}) the Cholesky decomposition to allow for fast
updating with new values of lambda, alpha and a.wght.  The "symbolic"
first step in the sparse Cholesky decomposition can also be used to
compute the sparse Cholesky decomposition for a different positive
definite matrix that has the same pattern of zeroes. This option is
useful in computing the likelihood under different covariance
parameters. For the LKrig covariance the sparsity pattern will be the
same if \code{NC}, \code{level}, \code{overlap} and the data locations
\code{x} are kept the same.  The returned component \code{LKinfo} has
class LKinfo and is a list with the information that describes the
layout of the multiresolution basis functions and the covariance
parameters for the GMRF.  (See help(LKinfo) and also \code{LK.basis}
as an example.)

\strong{\code{predict.LKrig, predictSE.LKrig: }} A vector of
predictions or standard errors.

\strong{\code{predictSurface.LKrig}: } A list in image format
(i.e. having components \code{x,y,z}) of the surface evaluated on a
regular grid. This surface can then be plotted using several different
R base package and fields functions e.g. \code{image},
\code{image.plot}\code{contour}, \code{persp}, \code{drape.plot}. The
\code{surface} method just calls this function and then a combination
of the image and contour plotting functions.



}
\author{
Doug Nychka
}

\seealso{ LatticeKrig, LKrig.sim.conditional,
mKrig, Krig, fastTps, Wendland, LKrig.coef, Lkrig.lnPlike, LKrig.MRF.precision,
 LKrig.precision
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


# breaking down the LatticeKrig function into several steps. 
# also use a different covariance model that has fewer absis fucntions
# (to make the example run more quickly)
 
  LKinfo<- LKrig.setup( x, nlevel=1, alpha=1, NC=15, a.wght=5,
                        lambda=.01)
# maximize likelihood over lambda see help( LKrig.MLE) for details
# this assumes the value of 5 for a.wght.  In many cases the fit is not
# very sensitive to the range parameter such as a.wght in this case --
# but very sensitive to lambda when varied on a log scale.

  MLE.fit<- LKrig.MLE(x,y, LKinfo=LKinfo)
  MLE.fit$summary # summary of optimization over lambda.
# fit using MLE for lambda MLE function has added MLE value of lambda to
# the LKinfo object.
  obj<- LKrig( x,y, LKinfo=MLE.fit$LKinfo.MLE)  
  print( obj)  

# find prediction standard errors at locations based on fixing covariance
# at MLE's
  out.se<- predictSE( obj, xnew= x)
# one could evalute the SE on a grid to get the surface of predicted SE's 
# for large grids it is better to use LKrig.sim.conditional to estimate
#  the variances by Monte Carlo

##########################################################################
# Use multiresolution model that approximates an exponential covariance
# Note that a.wght realted to a range/scale parameter
# is specified at a (seemingly) arbitrary value. 
##########################################################################
  
  LKinfo<- LKrig.setup( x, NC=6, nu=1, nlevel=3, a.wght= 5)
# take a look at the implied covariance function solid=along x
#  and  dashed=along y 
  check.cov<- LKrig.cov.plot( LKinfo)
  matplot( check.cov$d, check.cov$cov, type="l", lty=c(1,2))  

############################################################################
# Search over lambda to find MLE for ozone data with approximate exponential
# covariance
###########################################################################
\dontrun{
  
  LKinfo.temp<-  LKrig.setup( x, NC=6, nu=1,  nlevel=3, a.wght= 5)
# starting value for lambda optimzation
  LKinfo.temp$lambda<- 1 
  MLE.search<- LKrig.MLE(x,y, LKinfo=LKinfo.temp)
# this function returns an LKinfo object with the MLE for lambda included.
  MLE.ozone.fit<- LKrig( x,y,  LKinfo= MLE.search$LKinfo.MLE)
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
# single level with large range parameter -- similar to a thin plate spline
#  lambda specified 

# fit with elevation
  obj.CO.elev<- LKrig( x.CO,y.CO,Z=Z.CO, nlevel=1, NC=15, alpha=1, lambda=.005,
                          a.wght=4.1)
# BTW the coefficient for the linear term for elevation  is obj.CO$d.coef[4]
# fitted surface without the elevation term
\dontrun{
   LKinfo<- LKrig.setup( x.CO, nlevel=1, NC=20,alpha=1, a.wght=4.1, lambda=1.0)
  # lambda is the starting vlaue for MLE optimization
  CO.MLE<- LKrig.MLE( x.CO,y.CO,Z=Z.CO, LKinfo=LKinfo)
  obj.CO.elev<- LKrig( x.CO,y.CO,Z=Z.CO, LKinfo= CO.MLE$LKinfo.MLE)
  CO.surface2<- predictSurface( obj.CO.elev, drop.Z=TRUE, nx=50, ny=50)
# pull off CO elevations at same locations on grid as the surface
  data( RMelevation) # a superset of elevations at 4km resolution
  elev.surface<- interp.surface.grid( RMelevation, CO.surface2)
# linear adjustment for elevation 
  CO.surface.elev<- elev.surface$z*obj.CO.elev$d.coef[4]
  CO.full<- CO.surface2$z + CO.surface.elev
# for comparison fit without elevation as a linear covariate:
  CO.MLE2<- LKrig.MLE( x.CO,y.CO, LKinfo=LKinfo)
  obj.CO<- LKrig( x.CO,y.CO, LKinfo= CO.MLE2$LKinfo.MLE)
# surface estimate
  CO.surface<- predictSurface( obj.CO, nx=50, ny=50)
  set.panel( 2,1)
  coltab<- topo.colors(256)
  zr<- range( c( CO.full), na.rm=TRUE)
  image.plot( CO.surface,  col=coltab, zlim =zr)
    US( add=TRUE,lwd=2)
    title( "MAM min temperatures without elevation")
  image.plot( CO.surface$x, CO.surface$y, CO.full, col=coltab, zlim=zr)
    title( "Including elevation")
    US( add=TRUE,lwd=2)
}

########################################################################
# for a more elaborate search over  a.wght, alpha and lambda to find
# joint MLEs see help(LKrig.MLE)
########################################################################

########################################################################
# A bigger problem: 26K observations and 4.6K basis functions
# fitting takes about 15 seconds on a laptop for a fixed covariance
#  LKrig.MLE to find the MLE (not included) for lambda takes abou
#  8 minutes
#######################################################################
\dontrun{
  data(CO2)
  obj1<- LKrig( CO2$lon.lat,CO2$y,NC=100,nlevel=1, lambda=.088,
                       a.wght=5, alpha=1)
# 4600 basis functions 100X46 lattice  (number of basis functions
# reduced in y direction because of a rectangular domain
  obj1$trA.est # about 1040 effective degrees of freedom 
#
  glist<- list( x= seq( -180,180,,200),y=seq( -80,80,,100) )
  xg<-  make.surface.grid(glist)
  fhat<- predict( obj1, xg)
  fhat <- matrix( fhat,200,100) # convert to image
#Plot data and gap-filled estimate
  set.panel(2,1)
  quilt.plot(CO2$lon.lat,CO2$y,zlim=c(373,381))
  title("Simulated CO2 satellite observations")
  world(add=TRUE,col="magenta")
  image.plot( glist$x, glist$y, fhat,zlim=c(373,381))
  world( add=TRUE, col="magenta")
  title("Gap-filled global predictions")
}  
## same example but use a periodic lattice in longitude and 
## periodic basis functions. This is the "cylinder" model.
\dontrun{
  data(CO2)
  LKinfo.cyl<- LKrig.setup( cbind( c( -180,180), range( CO2$lon.lat[,2])),
                            NC=25,nlevel=1, lambda=.1,
                             a.wght=5, alpha=1,  distance.type="cylinder")
# 1127 basis functions from a 49X23 lattice
  search.CO2<- LKrig.MLE( CO2$lon.lat,CO2$y, LKinfo=LKinfo.cyl)
# MLE search over lambda
  obj2<- LKrig( CO2$lon.lat,CO2$y, LKinfo=search.CO2$LKinfo.MLE)
  surface( obj2)
  world( add=TRUE)
}  

set.panel()
#########################################################################
#  Comparing LKrig to ordinary Kriging
########################################################################

# Here is an illustration of using the fields function mKrig with the
# LKrig covariance to reproduce the computations of LKrig. The
# difference is that mKrig can not take advantage of any sparsity in
# the precision matrix because its inverse, the covariance matrix, is
# not sparse.  This example reinforces the concept that LKrig finds the
# the standard geostatistical estimate but just uses a particular
# covariance function defined via basis functions and the precision
# matrix.
# Load ozone data set (AGAIN)
\dontrun{
  data(ozone2)  
  x<-ozone2$lon.lat
  y<- ozone2$y[16,]
# Find location that are not 'NA'.
# (LKrig is not set up to handle missing observations.)
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
  a.wght<- 5
  lambda <-  1.5
  obj1<- LKrig( x,y,NC=16,nlevel=1, alpha=1,  lambda=lambda, a.wght=5,
                NtrA=20,iseed=122)
 
# in both calls iseed is set the same so we can compare 
# Monte Carlo estimates of effective degrees of freedom
  obj1$trA.est
# The covariance "parameters" are all in the list LKinfo
# to create this special list outside of a call to LKrig use
  LKinfo.test <- LKrig.setup( x, NC=16, nlevel=1, alpha=1.0,  a.wght=5)

# this call to mKrig should be identical to the LKrig results
# because it uses the LKrig.cov covariance with all the right parameters.
  obj2<- mKrig( x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
                      cov.args=list( LKinfo=LKinfo.test), NtrA=20, iseed=122)
# compare the two results this is also an
# example of how tests are automated in fields
# set flag to have tests print results
  test.for.zero.flag<- TRUE
  test.for.zero( obj1$fitted.values, obj2$fitted.values,
                  tag="comparing predicted values LKrig and mKrig")
# compare standard errors. 
  se1<- predictSE.LKrig( obj1)
  se2<- predictSE.mKrig(obj2)
  test.for.zero( se1, se2,
                  tag="comparing standard errors for LKrig and mKrig")
}
########################################################################
#  Unnormalized marginal variance for a 6X6 basis on [-1,1]X[-1,1]
#  This is an example of why normalization seems important. 
########################################################################

\dontrun{
# specify covariance without normalization note all we need is the
#corners of domains to setup the info list. 
  LKinfo<- LKrig.setup(cbind( c( -1,1), c(-1,1)), NC=6, nlevel=1,
                        a.wght=4.5,alpha=1, normalize= FALSE)  
# 80X80 grid of points 
  xg<- make.surface.grid( list(x=seq(-1,1,,80), y= seq(-1,1,,80))) 
  look<- LKrig.cov( xg, LKinfo=LKinfo,marginal =TRUE) 
# surface plot of the marginal variances of the process. 
  image.plot( as.surface(xg, look)) 
# basis function centers from the first (and only) level
  xp<- make.surface.grid( LKinfo$grid[[1]])
  points(xp)
}

}

\keyword{spatial}

