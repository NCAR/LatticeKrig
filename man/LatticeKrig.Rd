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
% # Modified  NL Jan -> DWN Jan 2014

\name{LatticeKrig}
\alias{LatticeKrig}
\alias{print.LatticeKrig}
\title{User-friendly spatial prediction and inference using a compactly supported
multi-resolution basis and a lattice model for the basis coefficients.
}
\description{ This is a simple and  high level function  to fit the LatticeKrig spatial
model to a data set. 
In is somplest form for 2-d spatial data:
\code{obj<-LatticeKrig(x,y)} will fit a LatticeKrig type model to 2d locations \code{x} and observation vector \code{y}. Several (actually many!) default choices are made for the multi-resolution spatial
covariance 
in this top level function. It uses the defaults that are based on 
a "thin-plate spline" like model for the spatial estimator  
 and also uses \code{LKrigFindLambda} to estimate  some covariance parameters (sill and nugget variances)
through likelihood maximization (i.e. estimates the measurement and process
variances.) For the simplest "black box" use only the observations and their 2-d
locations need to be specified. But please see the caveats below in the Details
section.

Despite the simple syntax,
LatticeKrig still takes advantage of the multi-resolution features of the basic 
\code{LKrig} function and any \code{LKrig} parameter can be passed through the
function call. See the example below for varying the range parameter.
Also, see \code{\link{LKinfo}} for documentation on the complete object that
describes the LatticeKrig model and \code{\link{LKGeometry}} for documentation on
extending or adding other spatial models to this package. 

The returned value from this function can be used subsequently for prediction, 
conditional simulation, and other parts of the spatial analysis. See 
\code{\link{predict.LKrig}} and \code{\link{LKrig.sim.conditional}} 
}


\usage{
LatticeKrig(x, y, Z = NULL, nlevel = 3, findAwght = FALSE, LKinfo = NULL,
 X=NULL, U=NULL, na.rm =
                 TRUE, tol = 0.005, verbose = FALSE, ...)
\method{print}{LatticeKrig}( x, digits=4, ...)
}

\arguments{
 \item{x}{Spatial locations of observations.  For the \code{LatticeKrig} function this should be a matrix
   where the columns index the spatial dimensions and rows index the observations.
For example for 100 2-d locations, \code{x} would be a  100X2 matrix. 

Or for the function \code{print.LatticeKrig} \code{x} is the returned object from the
    \code{LatticeKrig} function. }

  \item{y}{Spatial observations. No missing values are allowed.}
 
  \item{Z}{Linear covariates to be included in fixed part of the model
  that are distinct from the default first order polynomial in
  \code{x} (i.e. the spatial drift).}
  
     \item{X}{For linear inverse problems the matrix that maps coefficients of the basis to the   
  	 predicted values of observations. X must be in spam format. To convert from spind or dense
  	 format  to spam format see \code{help(spam)} as an alternative  
  	 \code{\link{spind2spam}}.}
  	 
  	 \item{U}{For linear inverse problems the matrix that maps coefficients of the   fixed part of the model to the predicted
	values of observations.}
  
   \item{nlevel}{Number of levels for the multi-resolution basis. Each
   level increases the number of basis functions by roughly a factor of
    4.}
    \item{findAwght}{If FALSE the default a.wght parameter (related to correlation range) is
    set to mimic a thin plate spline. If TRUE this parameter and hence the range is
    estimated my maximum likelihood. }
  \item{LKinfo}{An optional list giving the full specification of the
   covariance. If this is missing  it will be created internally and
   returned. If passed this parametrization will be used
   except lambda will be reestimated by maximum likelihood.}

  \item{na.rm}{If TRUE NA's are removed from \code{y} and \code{x} is subsetted.}

 \item{tol}{Tolerance for the log likelihood used to judge convergence.}

 \item{verbose}{ If TRUE print out intermediate results. }

  \item{\dots}{Additional arguments to pass to LKrig. See help(LKrig) for a complete list. But for convenience we note that  if you get some pesky memory warnings from spam
  you can set the storage higher by adding the argument \code{choleskyMemory}. For example to bump up to 2E6 include:
  \code{choleskyMemory=list(nnzR= 2E6)}.
  }
  
  \item{digits}{Number of significant digits in printed summary.}
}
\details{
Keep in mind that overall LatticeKrig is just a specific type of
spatial estimate that is designed to handle larger size data sets.
It focuses on a specific form of covariance function, but the estimator is
still the Kriging/Multivariate Gaussian Conditional Expectation/BLUE that is
standard in this field. 

The simplest model fit is:
 \deqn{Y_k = p(x_k)  + h(x_k) + e_k }
 \eqn{Y_k} is the \eqn{k^{th}} observation at location \eqn{x_k} with measurment error \eqn{e_k}.
 Here \eqn{p(x)} is a low order polynomial of degree m-1 with the default \code{m==2}. \eqn{h(x)} is a mean zero Gaussian process with the representation:
 \deqn{ h(x)= \sum_{l=1}^L \sum_{j=1}^{m(j)} \phi_{m,l}(x) c_{m,l} }
where \eqn{\phi} are multiresolution basis functions and the coefficients
have mean zero and spatial dependence specified by a Markov random field. Keep in mind that this unusual form still implies a specific covariance function for \eqn{h(x)}. In fact one can use the Krig or mKrig function from fields
to reproduce the LatticeKrig estimate for smaller data sets and check computations. 
(See \code{\link{LKrig.cov}} with examples ).  Details on the basis functions and the markov random field are given in the \code{\link{LKrig}}
help function.  
Throughout this package we assume the standard deivation of \eqn{e_k}
is \code{sigma}
and the marginal variance of  \eqn{h(x)} is \code{rho}. An important derived
parameter of the spatial estimate is \code{lambda = sigma^2/ rho} the noise to 
signal ratio. \code{sigma} and \code{rho} are estimated my restricted maximum
likelihood in \code{LatticeKrig}. 
 
This top level function is built on the more basic function \code{\link{LKrig}}
supports a very flexible
covariance.  \code{LKrig} depends on the parameters \code{nlevel}, \code{a.wght} and
\code{alpha}  specifying all these relevant parameters may be
discouraging to a new (or impatient!) user.  Thus LatticeKrig is a "wrapper" that
generates some simplified, default model choices to call the more general
function \code{LKrig}. It is useful for users not fully familiar with
the LatticeKrig methodology or those that wish to try a default
approach to get a quick look at a spatial analysis. You always go back and add some specific non default choices to the LatticeKrig call (e.g. changing \code{a.wght}). 
For the 2-dimensional  case the  default values
are set to give about 4 times as many basis functions as observations, use 5 extra
lattice points on the edge to minimize boundary effects, 
and to use four levels of multi-resolution. 
An important default is that a linear spatial drift is included in the model so
the model will relax to a linear prediction based on the x values in the absence of
a spatial component. In other words, the model includes by default a fixed part
that is linear in x.
The spatial correlation
range is nearly stationary and set large to mimic a thin-plate
spline. The smoothness mimics the Whittle covariance function (
smoothness = 1 for the Matern).  (See \link{LKrig.cov.plot} to get a plot of
the implied covariance function.) LatticeKrig also provides maximum likelihood 
estimates
of the measurement error standard deviation ("sigma") and process variance parameter ("rho") that are perhaps the 
parameters that most effect the shape of the estimated spatial field. The ubiquitous
parameter lambda throughout LatticeKrig is just the reparametrization lambda == sigma^2 / rho.

This top level function is pitched with all the caveats that
statistical model assumptions should always be checked and applying
generic methods to a specific problems without checking the
appropriateness can lead to misleading results. So plot your data and
try several models. Details on the full
computations can be found in the \code{LKrig} man page.
The \code{lambda = sigma2/ rho} parameter in \code{LKrig} is essential to the
Lattice Krig computation and an inappropriate value will result in
over or under fitting and incorrect interpolated values. The function
\code{LKrigFindLambda} is used within \code{LatticeKrig} to estimate a
lambda value from the data using maximum likelihood. 

One interesting feature of this package is the ability to handle spatial processes
on different geometries and the form is specified by the \code{LKGeometry} argument.
The current choices are:

\describe{
\item{\code{\link{LKRectangle}}}{A 2 dimensional Euclidean spatial domain. The default}
\item{\code{\link{LKInterval}}}{A 1 dimensional Euclidean spatial domain.}
\item{\code{\link{LKBox}}}{A 3 dimensional Euclidean spatial domain.}
\item{\code{\link{LKRing}}}{A 2 dimensional  spatial domain where the first coordinate is periodic (on [0,360]) and the second is Euclidean.  E.g. a slice around the equator and not having a large latitudal range. }
\item{\code{\link{LKCylinder}}}{A 3 dimension model where an additional  coordinate is added to the LKRing geometry. This is useful for representing a small section of the sphere where one also has a height component.}
\item{\code{\link{LKSphere}}}{A full 2-d spherical geometry. Coordiantes are given in longitude, latitude but the distances and any structures are on the sphere. 
}
}


One important feature of this package is that the different geometries all use the
same computation engine LKrig, following the same computational algorithm.  The
differences in setting up the problem and in evaluating the function are
implemented as S3 methods. The details of this strategy are described in 
\code{\link{LKGeometry}} and allow the user to add new geometries.  

This function also supports a model where the observations are simply
expressed as linear combinations of the basis function coefficients. Equivalently
this is a model where the observed data can be expressed as linear functinonals
applied to the polynmial term and the spatial process. 
Typically these 
linear maps represent observing integrals or weighted combinations of the 
fields and are important for data that is aggregated over by space.
See help(LKrig) for an example of how this model is set up at the end of the
Examples section.
}


\value{
The main call inside \code{LatticeKrig} is to \code{LKrig} and so a
\code{LKrig} object is returned. Thus all of the functionality that
comes with \code{LKrig} objects such as \code{predict},
\code{summary}, \code{predictSurface}, etc. remain the same as
described in \code{LKrig}.  Also, see the components \code{residuals}
and \code{fitted.values} in the returned object for these parts of the
fitted spatial model. The component \code{LKinfo} has all the details
that specify the basis functions and co variance model.
The component \code{MLE} gives details of the likelihood evaluations to estimate 
the sigma and rho parameters. 
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

# thin plate spline-like model with the lambda parameter estimated by
# maximum likelihood. Default choices are made for a.wght, nlevel, NC
# and alpha.

  obj<- LatticeKrig( x, y)
\dontrun{
# summary of fit and a plot of fitted surface
  print( obj)
  surface( obj )
  US(add=TRUE)
  points(x)
# prediction standard errors
  out.se<- predictSE( obj, xnew= x)
# predict at observations:
  out.fhat<- predict( obj, xnew= x)
# coveniently predict on a 100X100 grid for plotting
 out.surf<- predictSurface( obj, nx=100, ny=100)
# image.plot( out.surf) 
}
#
# In this example lon/lat are treated as just Euclidean coordinates 
# a quick adjustment to account for the difference in physical distance in N-S verses E_W
# is to just scale the longitude degrees to be comparable to degrees in latitude
# at least in the middle of the domain. The assumption is that for small spatial
# domains this approximation will not be bad for the coordinates at the edges too.
# You accomplish this by adding a scaling, V matrix:
\dontrun{
  meanLat<- mean( x[,2])*pi/180
  Vlonlat <- diag(  c( 1/cos(meanLat), 1) )
  obj1<- LatticeKrig( x, y, V = Vlonlat )
}

\dontrun{
# Refit using with just one level of  basis functions
# on a 20X20 grid witin the spatial domain ( so about 400) 
# actually number is 720 ( see obj1b$LKinfo) due adding edge nodes and also 
# aspect ratio of spatial domain 
# and find the a.wght parameter along with nugget and process variances. 
# this takes a while partly because LatticeKrig model is no optimized for small data sets!
  obj1b<- LatticeKrig( x, y, nlevel=1, NC=20, findAwght=TRUE)
# rudimentary look at how likelihood was optimized
#log lambda and omega =  log(a.wght-4)/2 are useful parametrizations ...
  quilt.plot( obj1b$MLE$lnLike.eval[,c("logLambda","omega")],
       obj1b$MLE$lnLike.eval[,"lnProfileLike.FULL"], 
       xlab="loglamda", ylab="omega",
       zlim =c(-640,-612))
  points( obj1b$MLE$lnLike.eval[,c("logLambda","omega")],cex=.25)
      
}
\dontrun{
# Refit using the tensor product type of basis functions
# (default is "Radial"). An example how an additional argument that is 
# passed to the LKrigsetup function to create the LKinfo object.
  obj2<- LatticeKrig( x, y, BasisType="Tensor")
}

#
# A 1-d example with 3 levels of basis functions
# See LKrig for an explanation if nlevel, NC,  alpha and a.wght 
# covariance parameters.


\dontrun{
 x<- matrix(rat.diet$t)
 y<- rat.diet$trt
 fitObj<- LatticeKrig( x, y)
# NOTE lots of defaults are set for the model! See print( fitObj)
 plot( x,y)
 xg<- matrix(seq( 0,105,,100))
 lines( xg, predict(fitObj, xg) )
}

\dontrun{
#  a 3D example
set.seed( 123)
N<- 1000
x<-  matrix( runif(3* N,-1,1), ncol=3, nrow=N)
y<-   10*exp( -rdist( x, rbind( c(.5,.5,.6) ) )/.5)

# NOTE setting of memory size for Cholesky. This avoids some warnings and
# extra computation by the spam package
LKinfo<- LKrigSetup( x,  nlevel=1,  a.wght= 6.01, NC=6, NC.buffer=2,
                    LKGeometry="LKBox", normalize=FALSE, mean.neighbor=200,
                    choleskyMemory=list(nnzR= 2E6) )                                      
out1<- LatticeKrig( x,y, LKinfo=LKinfo)

glist<- list( x1=seq( -1,1,,30), x2=seq( -1,1,,30), x3 = 0)
xgrid<- make.surface.grid( glist)

yhat<- predict( out1, xgrid)
# compare yhat to true function created above
image.plot( as.surface( glist, yhat))

}
#
###########################################################################
# Including a covariate (linear fixed part in spatial model)
########################################################################## 
\dontrun{
  data(COmonthlyMet)

  obj  <- LatticeKrig(CO.loc,  CO.tmin.MAM.climate, Z=CO.elev)
  obj2 <- LatticeKrig(CO.loc, CO.tmin.MAM.climate)

# compare with and without linear covariates
  set.panel(1,2)
  surface(obj)
  US(add=TRUE)
  title("With Elevation Covariate")

  surface(obj2)
  US(add=TRUE)
  title("Without Elevation Covariate")

}
\dontrun{
 data(COmonthlyMet)
# Examining a few different "range" parameters
a.wghtGrid<-  4  +  c(.05, .1, .5, 1, 2, 4)^2

#NOTE smallest is "spline like" the largest is essentially independent
# coefficients at each level.  In this case the "independent" end is
# favored but the eff df. of the surface is very similar across models
# indicating about the same separate of the estimates into spatial
# signal and noise
#
for( k in 1:5 ){
obj  <- LatticeKrig(CO.loc,  CO.tmin.MAM.climate, Z=CO.elev, 
                      a.wght=a.wghtGrid[k])
cat( "a.wght:", a.wghtGrid[k], "ln Profile Like:",
            obj$lnProfileLike, "Eff df:", obj$trA.est, fill=TRUE)
}

# MLE
obj0  <- LatticeKrig(CO.loc,  CO.tmin.MAM.climate, Z=CO.elev, 
                     findAwght=TRUE)
print(obj0$MLE$summary)
}

#########################################################################
# Reproducing some of the analysis for the example in the
# JCGS LatticeKrig paper.
#########################################################################

#### Here is an example of dealing with approximate spherical geometry.
\dontrun{
data(NorthAmericanRainfall)
library(mapproj)
x<- cbind(NorthAmericanRainfall$longitude, NorthAmericanRainfall$latitude)
y<- NorthAmericanRainfall$precip
log.y<- log(y)
elev<- NorthAmericanRainfall$elevation
# this is a simple projection as part of this and handled by the mapproj package
x.s<- mapproject( x[,1], x[,2], projection="stereographic")
x.s<- cbind( x.s$x, x.s$y)

# an alternative is to transform coordinates using another projection,
# e.g. a Lambert conformal projection
# with the project function from the rgdal package
# library( rgdal)
# x.s<- project(x,"+proj=lcc +lat_1=22 +lat_2=58 +lon_0=-93 +ellps=WGS84")
# this package has the advantage that the inverse projection is also 
# included ( inv=TRUE) so it is easy to evaluate the surface back on a Mercator grid.
             
obj0<- LatticeKrig(x.s, log.y, Z=elev )

fitSurface<- predictSurface( obj0, drop.Z=TRUE)
fitSurface$z<-  exp(fitSurface$z)/100
colorTable<- designer.colors( 256, c("red4", "orange", "yellow","green1", "green4", "blue"))
image.plot( fitSurface, col=colorTable)
map( "world", add=TRUE, col="grey30", lwd=3, proj="") 

}

}
\keyword{spatial}


