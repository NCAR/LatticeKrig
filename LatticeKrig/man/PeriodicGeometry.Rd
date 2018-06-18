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

\name{Spatial models for data on spherical regions. } 
\alias{LKRing}
\alias{LKCylinder}
\alias{LKSphere}
\title{Geometries to represent 
	2-d and 3-d spherical data. 
}
\description{These models include  a completely sphericial geometry based on a nearly regular node layout on the sphere. Also a simpler and perhaps more efficient 
  versions are implemented where the first coordinate is 
	periodic in the interval [0,360] and the remaining coordinates are regular (Euclidean).   
	This might be used to approximate a section of spherical data that excludes the 
	polar caps. These approximations are useful because one can take advantage of
	faster methods based on rectangular grids rather the more complex grids on a
	sphere. The disadvantage is that the mapping from these coordinates to the
	sphere is distorted as one gets close to the poles. 
  }
\section{Details}{
		These geometries are specified with the \code{LKGeometry} argument either in 
		LKrigSetup or LatticeKrig. The first coordinate is longitude either from [0,360]
		or [-180,180] and the second is lattiude [-90, 90].
		
		They each have the four specific methods: \code{LKrigLatticeCenters},
		\code{LKrigSAR}, \code{LKrigSetupLattice}, \code{setDefaultsLKinfo}
		and the source code is consolidated in the source files ModelRing.R and 
		ModelCylinder.R in the R subdirectory of this package source. For the spherical grid the
		a.wght is handled a but differently please read the note below.
		
		\bold{\code{LKSphere}} Recall that the core of the LatticeKrig model is the placement of
		the basis function centers and how the spatial autoregression is constructed from these
		node
		locations. Here the centers are generated according to a multiresolution based on the
		triangles from an Icsohedron. \code{\link{IcosohedronGrid}} creates a grid by taking the
		first
		level as the 12 points from a regular icsohedron. The subsequent levels generate a finer
		set
		of points by subdividing each triangular face into 4 new triangles. The three new mid
		points
		from the subdivision are added to the previous nodes to give the new level of resolution.
		The
		triangles tend to roughly  equilateral and so the nodes will tend to be roughly equally
		spaced
		. 
		\strong{NOTE:} Because the distances between nodes are not perfectly equally spaced  and nearest neighbors
		can be either 5 or 6 some adjustment is made to the
		weights when forming the SAR matrix. The net result is that it makes more sense to have the off
		diagonal rows sum to 1 and so the a.wght must be greater than 1.0 for a stationary model. 
		
		See the help file on 
		\code{link{IcosohedronFaces}} for code on visualzing these. The first 12 centers will
		have 5
		close neighbors and remaining centers will have 6. Currently the SAR weights are equal
		among
		all the neighbors and the \code{a.wght} parameter is adjusted by 
		\code{a.wght <- (a.wght - 6) + 5} for these locations. The basis functions have their 
		default as using great circle distance to determine the values between the center and
		other
		points. Because the distances between centers are not equal some adjustment is made to the
		\code{delta} parameter for the basis function support. The \code{NC} parameter in this
		geometry is the beginning resolution level. Currently the number of basis functions in
		each level are  
		\tabular{lrrrrrrrr}{
	  Level \tab   1 \tab     2   \tab   3  \tab    4  \tab    5  \tab
	  6   \tab   7   \tab   8 \cr
    Number Basis functions \tab  12 \tab    42 \tab   162 \tab   642 \tab 
    2562 \tab 10242  \tab 40962 \tab 163842 \cr
      }
So if \code{NC=2} and \code{nlevel=3} there will be a total of \code{42 + 162 + 642 = 846}
basis functions in the model. 
		
		\bold{\code{LKRing}}
		This model follows the Mercator projection for a sphere where 
		longitude and latitude are treated as Euclidean coordinates except that 
		longitude is periodic. So the actual coordinates represent the surface of 
		 cylinder which is one way of visualizing the Mercator projection. 
		To keep things simple the first coordinate is 
		essentially hardwired to be in the scale of degrees (sorry for all you fans of radians)
		and wrapping 
		0 to 360 ( or periodic in [-180,180]). It is important to scale the second coordinate in
		this
		geometry to be comparable in spatial scale to degrees (use the \code{V}
		argument in LKrigSetup). However, if the second coordinate can be interpreted as a
		latitude it
		is often reasonable to assume the spatial scales are the same in these two coordinates. 
		
		Note this geometry can also be used to represent an equatorial section of a spherical
		volume.
		Here the first coordinate is longitude but the second can be 
		interpreted as a radius from the sphere's center. 
		This is a case where care needs to taken to scale the radial component to make sense with
		the degrees in the first.  
		
		\bold{\code{LKCylinder}}
		This is just the three dimensional extension of LKRing
		with the first coordinate being periodic in (0,360) and the
		remaining two treated as Euclidean	
	
	The periodicity in the first coordinate is implemented in two places. 
	First in setting up the spatial autoregression (SAR) weights, the weights
	reflect the wrapping. Second in finding distances between coordinates the
	nodes in the lattice has the first coordinate tagged as being periodic. 
	Specifically in LKrigSetupLattice the gridList for each lattice has an attribute vector that
	indicates which coordinates are periodic. This information is used in the distance function
	LKrigDistance when called with arguments of a matrix and a gridList. 
	
	Why is this so complicated? This structure is designed around the fact that one 
	can find the pairwise distance matrix quickly between an arbitrary set of locations and a
	rectangular grid (a gridList object in this package).
	The grid points within a delta radius of an arbitrary point can be found by simple
	arithmetic
	and indexing.  Because these two geometries have regular 
	lattice spacings is it useful to exploit this. See \code{ \link{LKrigDistance}}
	for more details about the distance function. 
	
	Finally, we note that for just patches of the sphere one can use the usual
	LKRectangle geometry and change the distance function to either chordal or 
	great circle distance. This gives a different approach to dealing with the
	inherent curvature but will be awkward as the domain is close to the poles. 
			
	}
	
\author{
	Doug  Nychka and Zachary Thomas 
}
\examples{
# 	
# fit the CO2 satellite data with a fixed lambda
# (use a very small and unrealistic number of basis functions and levels so example
#  runs quickly)
 data(CO2)
# to look at raw data: quilt.plot(CO2$lon.lat, CO2$y, nx=288, ny=165)
# NOTE: NC=1 is related to how grid for sphere is set up.
  LKinfo0<- LKrigSetup( CO2$lon.lat, NC=1 ,nlevel=2,
                       a.wght=6.1, alpha=c(1,.25),
                       LKGeometry="LKSphere" )                                         
  # Use MLE for lambda found below 
  obj0<- LKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo0, lambda=0.01068)
  surface(obj0, nx=288, ny=165)
  world( add=TRUE)
\dontrun{
# estimate lambda ( should be 0.01068)
obj0B<-  LatticeKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo0)
# a more serious model this uses about 3300 basis functions
LKinfo0<- LKrigSetup( CO2$lon.lat, NC=3 ,nlevel=3,
                       a.wght=6.1, alpha=c(1, .5, .25),
                       LKGeometry="LKSphere" )
                       
obj0B<-  LatticeKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo0)
# takes about 1 minute on a macbook air
# setting findAwght = TRUE  takes about 8 minutes with 
# lambda = 1.737 and a.wght = 15.8 (independence among basis coefficients)
} 
# a more serious model 
\dontrun{
  LKinfo1<- LKrigSetup(CO2$lon.lat, NC=8 ,nlevel=1, lambda=.2,
                       a.wght=5, alpha=1, 
                       LKGeometry="LKRing" )                                         
  obj1<- LatticeKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo1)	
# take a look: 
surface( obj1)
world( add=TRUE) 
}
# compare to fitting without wrapping
#  LKinfo2<- LKrigSetup(CO2$lon.lat, NC=8 ,nlevel=1,
#                   lambda=.2, a.wght=5, alpha=1 )                                         
#  obj2<- LKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo2)	
# not periodic in longitude:
# surface(obj2)  

# a synthetic example and larger example
\dontrun{
 set.seed(124)
 N<- 1e4
 NC<- 3
  x0<- matrix( rnorm(3*N), ncol=3)
  x0<- x0/ sqrt( rowSums( x0^2))
  
  x<-  toSphere( x0 )
  
  # the true function for testing -- a bump at the direction alpha
  fun<- function(X){
    alpha<-  c( .4,.2,1)
    alpha<- alpha/ sqrt( sum( alpha^2))
    4*( 1 + c(( X)\%*\%alpha) )^2 
  }
  
  ytrue <- fun(x0)
  y<- ytrue + .05*rnorm( length(ytrue))
  xr=cbind( c(-180, 180), c(-90,90))
  # this defines about 3300 basis functions
  LKinfo1<- LKrigSetup( xr,
                        NC=NC, LKGeometry="LKSphere", a.wght=1.01,
                        nlevel=3, alpha = c(1.0,.5,.25)^2,
                        choleskyMemory=list(nnzR= 20E6),
                        normalize=FALSE)
  out<- LatticeKrig( x,y, LKinfo=LKinfo1)                      
surface( out)                        
}

}

\seealso{
	\code{\link{LatticeKrig}},
  \code{\link{LKrigSetup}},
	\code{\link{LKrigSAR}},
	\code{\link{LKrigLatticeCenters}}
	}
	


\keyword{ spatial }
 
