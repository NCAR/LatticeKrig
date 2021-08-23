\name{directionCosines}
\alias{directionCosines}
\alias{toSphere}
\alias{projectionSphere}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Utility functions for spherical coordinate and projections. 
%%  ~~function to do ... ~~
}
\description{
Convert back and forth between  lon/lat and direction cosines and also project spherical coordinates on a local tangent plane. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
directionCosines(x)
toSphere( Grid)
projectionSphere(x1, x2) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A two column matrix of lon/lat coordinates in degrees.
%%     ~~Describe \code{x} here~~
}
\item{x1}{ A vector of direction cosines defining the intersection of the tangent plane
with the unit sphere}
\item{x2}{ A matrix of direction cosines that will be projected onto the plane defined by 
 x2.}
\item{Grid}{
A three column matrix of direction cosines 
%%     ~~Describe \code{x} here~~
}
}
\details{
The conversion functions are  based on straight forward definitions of spherical 
coordinates.  The projection function is done by two rotations, first around the Z axis and then around the Y axis.  
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\strong{directionCosines} A three column matrix of direction cosines.

\strong{toSphere} A two column matrix of longitudes and latitudes.

\strong{projectionSphere} A two column matrix of Euclidean coordinates on the tangent plane to \code{x1}. The convention is that the origin (0,0) is mapped to \code{x1} and the X- axis are points along the meridian through  \code{x1}. The Y axis are points on the great circle passing through  \code{x1} and perpendicular to the meridian. 

}

\author{
Doug Nychka
%%  ~~who you are~~
}

\examples{
# 
# icosahedron:
 phi = (1 + sqrt(5))/2
  V = matrix(c(0, 1, phi, 0, 1, -phi, 0, -1, phi, 0, -1, -phi, 
        1, phi, 0, -1, phi, 0, 1, -phi, 0, -1, -phi, 0, phi, 
        0, 1, -phi, 0, 1, phi, 0, -1, -phi, 0, -1), 12, 3, byrow = TRUE)

# check : library( rgl); plot3d( V, size=10, col="red4" )
# as lon lat:
V2<- toSphere( V)
plot( V2)

# lon lat grid
 lGrid<- make.surface.grid(  list(x= seq( -10,10,, 10), y= seq( -20,20,,10)) )
 
 dGrid<- directionCosines( lGrid)
 pairs( dGrid)
 # also try:   library( rgl); plot3d(  dGrid)
 }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

