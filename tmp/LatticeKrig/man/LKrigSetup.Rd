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

\name{LKrigSetup}
\alias{LKrigSetup}
\alias{LKinfoUpdate}
\alias{LatticeKrigEasyDefaults}
\alias{LKinfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create or update the LatticeKrig model object  (LKinfo) for spatial fitting. 
}
\description{
This function combines some input arguments with defaults for other
to create a list 
describing the LatticeKrig spatial model. A key to specifying the LatticeKrig spatial model is specifying the geometry, e.g. LKRectangle for a 2-d rectangular domain. Each geometry has some parameters that control the basic model setup and these are included through the \dots arguments of  this function. See the help for this argument below for some examples.
}
\usage{
LKrigSetup(x = NULL, nlevel = NULL, alpha = NA, alphaObject =
                 NULL, nu = NULL, a.wght = NA, a.wghtObject = NULL, NC
                 = NULL, NC.buffer = NULL, normalize = TRUE, lambda =
                 NA, sigma = NA, rho = NA, rho.object = NULL,
                 latticeInfo = NULL, basisInfo = NULL, LKGeometry =
                 "LKRectangle", distance.type = "Euclidean",
                 BasisFunction = "WendlandFunction", overlap = 2.5, V =
                 NULL, BasisType = "Radial", fixedFunction =
                 "LKrigDefaultFixedFunction", fixedFunctionArgs =
                 list(m = 2), collapseFixedEffect = FALSE, max.points =
                 NULL, mean.neighbor = 50, choleskyMemory = NULL,
                 verbose = FALSE, noCheck = FALSE, returnCall = FALSE,
                 dense = FALSE, ...)
                 
LatticeKrigEasyDefaults(argList, nlevel, x)

LKinfoUpdate(LKinfo, ... )
}

\arguments{	

\item{argList}{Argument supplied to the top level LatticeKrig function.}
	
\item{alpha}{A vector of length \code{nlevel} with the relative
variances for the different multi-resolution levels.  }
\item{alphaObject}{For non-stationary models an object to be used with the predict function to give the alpha values at the process locations.
Typically this is a list of "predict" objects. 
See \link{nonstationaryModels} for examples how to use this option. 
}

 \item{a.wght}{This parameter controls the correlation range in the SAR
 model. In most cases, a scalar value and for the 2-d model by default
 greater than 4. 
   If a vector this can specify an anisotropic set of weights. To
   specify a.wght parameters that are different for each level they
   should be in the form of a list 
   of length \code{nlevel}. E.g.  \code{a.wght = list( 4.5, 5, 10)}
   will specify different a.wghts for 3 different levels.   The setup
   function will
    check that the values are in a valid range for a geometry and the
    length of the list agrees with the number of levels. 
    
  	The details of how this is connected to the covariance function
  	varies
  	based on the geometry. However, qualitatively this is related to a
  	range parameter. For the LKRectangle geometry and a stationary
  	model,  at level k the center point has
  	weight 1 with the 4 nearest neighbors given weight
  	-1/a.wght[k]. In this case a.wght must be greater than 4 for the
  	fields to be
  	stationary and following Lindgren and Rue the range parameter is
  	approximately \code{1/sqrt(a.wght-4)}. 
}

\item{a.wghtObject}{For non-stationary models an object to be used with the predict function to give the a.wght values at the lattice locations. See \link{nonstationaryModels} for examples how to used option. 
}

 \item{basisInfo}{ A list with extra components the object used to
  	describe the multi-resolution basis.  Usually this will not be
  	needed for standard models. }
  	
  \item{BasisType}{A character string indicating the type of basis function. Currently
  	this is either \code{"Radial"} or \code{"Tensor"}. }

\item{choleskyMemory}{A list that will be used in the spam call to 
do the Cholesky decomposition. See the \code{memory} argument in 
\code{chol.spam}. }

\item{collapseFixedEffect}{If FALSE the fixed part of the model is
found separately for each replicated data set. If TRUE the estimate is
polled across replicates.This is largely a modeling decision whether
variation among the replicate fields is due to the spatial component or
also include variation in the fixed effects across replicates -- guess
they are not really fixed then!}

\item{dense}{If FALSE sparse linear algebra is used for the
computations
. If TRUE the matrices are
made "dense" (zeroes are filled in) and the ordinary Lapack functions
are used for the linear algebra.
This option is primarily
for testing and timing sparse verses standard algorithms. }

 \item{distance.type}{A text string indicate type distance to use
  	between spatial locations when evaluating the basis
  	functions. Default is "Euclidean". Other choices when
  	locations are in degrees longitude and latitude are "Chordal"
  	and "GreatCircle".  The default radius is in miles (3963.34)
  	but can be reset using the "Radius" attribute to the text
  	string.

  	e.g.  \code{ dtype<- "Chordal"; attr( dtype, which="Radius")
  		  <- 6371 }

  	will set radius to kilometers. }
  	
\item{fixedFunction}{A text string that is the name of the function
	used to find the fixed part of the spatial model based on the
	locations. The default is a linear (m=2) polynomial in the
	spatial coordinates. See \code{LKrigDefaultFixedFunction} for more details.
	Set this to NULL if you do not want to
	include a fixed effect in the spatial model. }
	
\item{fixedFunctionArgs}{A list containing arguments to supply when
	evaluating the fixed function.}
	
    \item{lambda}{ The "noise to signal ratio" or also known as the
  	smoothing parameter it is the parameter \code{lambda =
  	sigma^2/rho}. If specified then sigma and rho typically are
  	estimated in LKrig by maximum likelihood. If lambda is not
  	specified then it is set as \code{lambda = sigma^2/ rho}. Note
  	that to evaluate the spatial process model, e.g. using the
  	function \code{LKrig.cov}, a value of lambda is not needed and
  	this argument can default to NA.  } \item{latticeInfo}{ Part
  	or all of the object used to describe the Markov random field
  	lattice.  In the standard cases this list is created in the
  	setup function and need not be specified. See
  	\code{LKrigSetupLattice} for details. Note that the contents
  	of this list is concatenated to any additional components
  	supplied by \code{LKrigSetupLattice}.}
  	 
  \item{LKGeometry}{ A text string that gives the names of the model
  	geometry.  The default is "LKrectangle" assuming the 
  	spatial domain is a rectangle.  Other common choices are
  	"LKInterval" (1 d problem ) and "LKBox" (3 d problem) }

\item{LKinfo}{A list that has class "LKinfo".}

  \item{mean.neighbor}{ The average number of nonzero points when each
        basis function is evaluated at a set of points points in the spatial
        domain.}
 
  \item{max.points}{This is a parameter for the nearest neighbor
  	distance function that sets the maximum array size for the
  	nonzero distances among points.  e.g. with 100 points each
  	with 20 nonzero neighbors max.points needs to be 2000 =
  	100*20. Specifically if the total number of nonzero values
  	when each basis function is evaluated at all the spatial
  	locations. The easier way to specify space is by using the
  	mean.neighbor argument.  } 
  	
\item{NC}{For regular grids of lattice points the maximum number of lattice grid points for a spatial coordinate and at the coarsest level of resolution. For a example, for a square region,
(and \code{LKGeometry = "LKRectangle"})
\code{NC=5} results in a 5X5  = 25  lattice points at the first level.
Note that the default is that lattice points are also taken to have the same spacing in every dimension. } 
\item{NC.buffer}{Number of extra lattice points added outside the spatial domain 
for regular grids of lattice points. This helps to reduce boundary effects from the
SAR model.  \code{NC.buffer=5} and \code{NC=5} for a square region will result in 
(5+ 2*5) X (5+ 2*5) = 225 lattice locations at the coarsest level of resolution.
Note that this number by default is fixed for finer resolutions and so does not
contribute as much to the total number of lattice points.}
 \item{nlevel}{Number of levels in
  	multi-resolution. Note that each subsequent level increases the
  	number of basis functions within the spatial domain size by a
  	factor of roughly 4. }
 
 \item{noCheck}{If FALSE do not make any checks on the consistency of
  	the different parts of the final LKinfo object. e.q. values of
  	a.wght within the range of a stationary Markov random field. }
 
  \item{normalize}{ If TRUE the basis functions will be normalized to
  	give a marginal variance of one for each level of
  	multi-resolution.  (Normalizing by levels makes it easier to
  	interpret the \code{alpha} weights.)  } 
  
  \item{nu}{ A
  	smoothness parameter that controls relative sizes of
  	\code{alpha}.  Currently this parameter only makes sense for
  	the 2D rectangular domain (LKRectangle) } 
  	
  \item{overlap}{
  	Controls the overlap among the radial basis functions and
  	should be in units of the lattice spacing. For the rectangular
  	case the default of 2.5 means that the support of the Wendland
  	basis function will overlap 2.5 lattice units in each
  	direction.  See \code{LKrig.basis} for how overlap adjusts
  	scaling in the basis function formula. }
  	
  \item{rho}{ A scalar, the sill or marginal variance of the process.
  	} 
  	
  \item{BasisFunction}{Text string giving the 1-d form for
  	basis function. }

  \item{rho.object}{ A prediction object to specify part of the
  	marginal variance for the process. Specifically the form is
  	VAR(g(x1))= rho*h(x1) Calling \code{predict(rho.object,x1)}
  	should return a vector with the values of h at the (arbitrary)
  	spatial locations in x1.  If omitted assumed to be the
  	constant one -- the usual case.  }
 	\item{returnCall}{If TRUE the call to LKrigSetup is also included as
 	      part of the LKinfo object. }
 
  \item{sigma}{ The measurement error standard deviation.  Also called
  	the nugget variance in geostatistics.  }
  	
  \item{V}{See entry in \code{LKrig}.}
 
  \item{verbose}{If TRUE print out intermediate information.}
  
  	\item{x}{Spatial locations that define the spatial domain for
prediction. This is only used to determine ranges of the grid for the basis
functions so, for example, for a rectangular domain only two points
are required that bound the rest of the data locations. E.g.  x=
rbind( c( 0,0), c(1,1)) will set the domain to be the unit square.  }

 	
  \item{\dots}{Specific arguments that will be included in the
  	\code{setupArgs} list and also passed to
  	\code{LKrigSetupLattice}. For \code{LKinfoUpdate} these specify the components
  	       of LKinfo to update. For example for a rectangular domain (the default) the argument \code{NC} is needed to specify the initial lattice size \code{NC.buffer} sets the number of extra points included in the lattice and defaults to 5.  \code{NC} for the LKSphere geometry  is the initial level of the
	       icosahedral grid.   }
}
\details{
	Many of the functions within LKrigSetup are overloaded to
	adapt to the \code{LKGeometry} class.  This makes it easy to add new
	geometries or other models to the LatticeKrig framework.  The
	required components of this object (see below) outline how
	the latticeKrig model is structured and what should be common
	features independent of the geometry. The key components are
	latticeInfo that contains the information used to generate the
	spatial autoregressive matrix on the lattice (see
	\link{LKrigSAR}) and the lattice centers (see
	\link{LKrigLatticeCenters}). The component basisInfo used to generate
	the radial basis function (see LKrig.basis)

	The function LKrigEasyDefaults is used in the top level 
	function LattticeKrig to make the logic of different default choices 
	more readable and reduces the clutter in this function. Its main 
	purpose is to find a reasonable choice for NC when this is not
	specified. 

	The function LKinfoUpdate is more of a utility used for
	clarity that allows one to update the LKinfo object with
	particular components without having to recreate the entire
	object. This function is used in the MLE search when just
	values of alpha or a.wght are being varied.
}

\value{
	An object with class "LKinfo" and also the additional class given by LKGeometry. The required components are:
\describe{
	\item{nlevel}{ Number of levels}
	\item{alpha}{alpha parameters as a list that has  nlevel components and possibly some attributes. }
	\item{a.wght}{a.wght parameters as a list that has  nlevel components and possibly some attributes. }
	\item{nu}{nu parameter}
	\item{normalize}{A logical indicating whether to normalize.}
	\item{lambda}{ Value of lambda.}
	\item{sigma}{Value of sigma.}
	\item{rho}{Value for rho.}
	\item{rho.object}{Value for rho.object.}
	\item{latticeInfo}{A list with specific multi-resolution lattice information}
  \item{setupArgs}{All arguments passed in the call and any in in \dots }
  \item{basisInfo}{A list with basis information.}
  \item{call}{The actual call used to create this object.}	 
}

}

\author{
	Doug Nychka
}

\examples{
  data(ozone2)
  # find the ranges of the  data, this is the same as passing
  # the entire set of observation locations and is more compact 
  sDomain<-apply( ozone2$lon.lat, 2,"range")
  LKinfo0<- LKrigSetup( sDomain, NC=10, nlevel=2, alpha=c(1,.5),
                       a.wght = 5)
  print( LKinfo0)
  
  #Gigantic buffer added note extra basis functions. 
  LKinfo<- LKrigSetup( sDomain, NC=10, NC.buffer= 15, nlevel=2, 
  alpha=c(1,.5),a.wght = 5)
  print( LKinfo)
  
  LKinfo2<- LKinfoUpdate( LKinfo,  a.wght=4.1, NC=12)
  LKinfo3<- LKrigSetup( sDomain, NC=12, nlevel=2, alpha=c(1,.5),
                        a.wght=4.1)
# LKinfo2 and LKinfo3 should be the same.                         
 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial}

