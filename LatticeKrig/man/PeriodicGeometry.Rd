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
\description{These models include  a completely spherical geometry based on a nearly regular node layout on the sphere. Also a simpler and perhaps more efficient 
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
		or [-180,180] and the second is latitude [-90, 90].
		
		They each have the four specific methods: \code{LKrigLatticeCenters},
		\code{LKrigSAR}, \code{LKrigSetupLattice}, \code{setDefaultsLKinfo}
		and the source code is consolidated in the source files ModelRing.R and 
		ModelCylinder.R in the R sub-directory of this package source. For the spherical grid the
		a.wght is handled a but differently please read the note below.
		
		\bold{\code{LKSphere}} Recall that the core of the LatticeKrig model is the placement of
		the basis function centers and how the spatial autoregression is constructed from these
		node
		locations. Here the centers are generated according to a multi-resolution based on the
		triangles from an icosahedron. \code{\link{IcosahedronGrid}} creates a grid by taking the
		first
		level as the 12 points from a regular icosahedron. The subsequent levels generate a finer
		set
		of points by subdividing each triangular face into 4 new triangles. The three new mid
		points
		from the subdivision are added to the previous nodes to give the new level of resolution.
		The
		triangles tend to roughly equilateral and so the nodes will tend to be roughly equally
		spaced. This regularly improves for the finer generations of triangles, however,
		the 12 original nodes from the icosahedron will have 5 nearest neighbors rather than 6. 
		
	  The setup argument \code{startingLevel} specifies the first level of the lattice and \code{nlevel} the
	  total number. 
	  
		\strong{NOTE:} Because the distances between nodes are not perfectly equally spaced
		and nearest neighbors
		can be either 5 or 6 some adjustment is made to the
		weights when forming the SAR matrix. The net result is that it makes more sense to have the off
		diagonal rows sum to 1 and so 
		\strong{ a.wght must be greater than 1.0 } for a stationary model.                
		
		See the help file on 
		\code{link{IcosahedronFaces}} for code on visualizing these. The first 12 centers will
		have 5
		close neighbors and remaining centers will have 6. Currently the SAR weights are
		roughly equal
		among
		all the neighbors but are adjusted so that a locally linear function will be in the
		"null" space of these weights. For each node the neighbors are projected onto the tangent
		plane to the sphere at this node location. Now consider a linear function on the coordinates
		in this tangent plane. The idea is to find  weights applied to the neighbors that will give
		a perfect linear prediction for the value at the node. The negatives of these values are
		used as the SAR weights. Specifically let \eqn{w_k} be the weights and \eqn{Y_k} the values of the
		field at the these locations, and \eqn{Y_0} the value at the node.  Then by construction
\deqn{ Y_0 - sum w_k Y_k  = 0 }
                for any field that is linear in the tangent plane to the node. Specifically in the function
		\code{LKrigSAR.LKSphere} the weights follow the code fragment
		
		\preformatted{
		   x1<- grid3d[J,]
       x0<- grid3d[I,]
       u<- projectionSphere( x0,x1) 
    # u are local 2 d coordinates on tangent plane to sphere at x0
    # x0 projects to (0,0)
      X<- cbind( rep( 1,nJ), u )
    # find weights to predict a linear function
    # at node from the nearest neighbors. 
      W<- c( (X)%*%(solve( t(X)%*%X, c( 1,0,0) )  ))
    }
		 \code{I} is the index of the node at lattice point \code{x0},
		 \code{J} is the index of the neighbors at lattice points \code{x1},
		 and \code{-W} the weights used in the SAR.
		 
		
		The basis functions have their 
		default as using great circle distance to determine the values between the center and
		other
		points. Because the distances between centers are not equal some adjustment is made to the
		\code{delta} parameter for the basis function support. Currently the number of basis functions in
		each level are  
		\tabular{lrrrrrrrr}{
	  Level \tab   1 \tab     2   \tab   3  \tab    4  \tab    5  \tab
	  6   \tab   7   \tab   8 \cr
    Number Basis functions \tab  12 \tab    42 \tab   162 \tab   642 \tab 
    2562 \tab 10242  \tab 40962 \tab 163842 \cr
      }
So if \code{startingLevel=2} and \code{nlevel=3} there will be a total of \code{42 + 162 + 642 = 846}
basis functions in the model if the spatial domain includes the entire sphere. 
		
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
# (but use a very small, unrealistic number of basis functions and levels so example
#  runs quickly)
 data(CO2)
# to look at raw data: quilt.plot(CO2$lon.lat, CO2$y, nx=288, ny=165)
# Use two levels starting at the second generation of lattice points from the triangulation
  LKinfo0<- LKrigSetup( CO2$lon.lat, startingLevel=1 ,nlevel=2,
                       a.wght=1.1, alpha=c(1,.25),
                       LKGeometry="LKSphere" )
# Take a look at Model summary
  print( LKinfo0)
  
# Use arbitrary  lambda  (.01) 
  obj0<- LKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo0, lambda=0.01)
\dontrun{
# Surface plot of estimate
  surface(obj0, nx=288, ny=165)
  world( add=TRUE)
}  
\dontrun{
data(CO2)
# estimate lambda ( should be around 0.003)
# NOTE: lambda values will tend to be sensitive to the model choice
  LKinfo0<- LKrigSetup( CO2$lon.lat, startingLevel=2 ,nlevel=2,
                       a.wght=1.1, alpha=c(1,.25),
                       LKGeometry="LKSphere") 
  obj0B<-  LatticeKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo0)
# a more serious model this uses about 3300 basis functions
LKinfo0<- LKrigSetup( CO2$lon.lat, startingLevel=3, ,nlevel=3,
                       a.wght=1.1, alpha=c(1, .5, .25),
                       LKGeometry="LKSphere" )
                       
obj0B<-  LatticeKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo0)
# takes about 1 minute on a Macbook air
# setting findAwght = TRUE  takes about 8 minutes with 
# lambda = 1.737 and a.wght = 15.8
} 
#####################################
# The ring geometry
#####################################
\dontrun{
  data(CO2)
  LKinfo1<- LKrigSetup(CO2$lon.lat, NC=8 ,nlevel=1, lambda=.2,
                       a.wght=5, alpha=1, 
                       LKGeometry="LKRing" )                                         
  obj1<- LatticeKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo1)	
# take a look: 
surface( obj1)
world( add=TRUE) 
}
# compare to fitting without wrapping:
\dontrun{
  LKinfo2<- LKrigSetup(CO2$lon.lat, NC=8 ,nlevel=1,
                   lambda=.2, a.wght=5, alpha=1 )                                         
  obj2<- LKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo2)	
 # NOTE: not periodic in longitude:
 surface(obj2)  
}

# a synthetic example and larger example
\dontrun{
 set.seed(124)
 N<- 1e4
  x0<- matrix( rnorm(3*N), ncol=3)
  x0<- x0/ sqrt( rowSums( x0^2))
  
  x<-  toSphere( x0 )
  
# the true function for testing -- a bump at the direction alpha
  fun<- function(X){
    alpha<-  c( .1,.1,1)
    alpha<- alpha/ sqrt( sum( alpha^2))
    4*( 1 + c(( X)\%*\%alpha) )^2 
  }
  
  ytrue <- fun(x0)
  y<- ytrue + .05*rnorm( length(ytrue))
# this defines about 3300 basis functions
  LKinfo1<- LKrigSetup( x,
                        startingLevel=3,
                        LKGeometry="LKSphere",
                        a.wght=1.01,
                        nlevel=3, alpha = c(1.0,.5,.25)^2,
                        choleskyMemory=list(nnzR= 20E6),
                        normalize=TRUE)
  out<- LatticeKrig( x,y, LKinfo=LKinfo1, lambda=.01)                      
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
 
