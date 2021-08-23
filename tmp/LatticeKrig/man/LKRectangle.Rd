\name{LKRectangle}
\alias{LKRectangle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary of the LKRectangle geometry for a standard two
dimensional spatial domain.
}
\description{
The basic LatticeKrig model is for a  2-d spatial
domain and is identified by the geometry class
LKRectangle. The lattices used in this case are equally spaced 
regular grids that are nested and the default distance function
is standard Euclidean distance. The number of lattice points in
the first level is NC inside the spatial domain and the subsequent
levels decrease the lattice spacing by factors of 2.  The spatial
autoregressive
coefficients, referred to as the a.wght parameter(s) have several levels
of detail depending on the
stationarity and isotropy. In the simplest case a.wght is the central
value, and should be greater than 4. Also the four nearest neighbors
are take to be -1. In the most complicated \code{a.wght} can include all 8
nearest 
neighbors can be specified differently at every lattice node. To be stable,
the sum of the 
\code{a.wght} parameters for a node should be greater than zero. 
}
\details{
Here is a simple and small three level example that sets up the LatticeKrig model
for spatial estimation and prediction. It assumes a spatial domain with
extent [-1,1] in the horizontal and [-1,1] in vertical and beginning 
with 4 grid points in lowest level and with a default of 5 extra points
outside the domain for boundary corrections. The assumed SAR model defaults
to a central value of 4.1 and the 4 nearest weights are -1. Defining the
extent of the spatial domain explicitly was done to simplify the example. Typically
one just defines the domain to the minimum and maximum x and y locations 
of the spatial data (and the locations can passed as the first argument 
in the example below instead of finding ranges.) Note that does not mean
that the spatial locations
fill out a rectangle and they do not need to regularly spaced. The \code{nu} parameter is a handy way to specify the relative weights given each level. For a larger numbers of levels this parameter is equivalent to the Matern smoothness parameter. 

\strong{Setting up the LKinfo object.}

\preformatted{
 sDomain<- cbind( c(-1,1), c( -1, 1))
LKinfo<- LKrigSetup(sDomain, nlevel=3, NC=4, a.wght=4.1, nu=1.0)
}

with the result 

\preformatted{
print(LKinfo)
Classes for this object are:  LKinfo LKRectangle
The second class usually will indicate the geometry
     e.g.  2-d rectangle is  LKRectangle
 
Ranges of locations in raw scale:
     [,1] [,2]
[1,]   -1   -1
[2,]    1    1
 
Number of levels: 3
delta scalings: 0.6666667 0.3333333 0.1666667
with an overlap parameter of  2.5
alpha:  0.7619048 0.1904762 0.04761905
based on smoothness nu =  1
a.wght:  4.1 4.1 4.1
 
Basis  type: Radial using  WendlandFunction  and 
Euclidean  distance.
Basis functions will be normalized
 
Total number of basis functions  1014
 Level Basis size      
     1        196 14 14
     2        289 17 17
     3        529 23 23
 
Lambda value:  NA
}

\strong{About the lattice}
The number of nodes define the number of basis
functions and at first may seem a bit mysterious. 
However at the first level we get 4 lattice points
with 5 extra boundary points added in the x
direction amounting to 14 total and similarly 14
in the y direction because it is exactly of the
same size. The default is the the parameter NC
determines the number of lattice points in the larger
dimension and the smaller dimension is divided according 
to the spacing from the longer dimension. For example if 
the y extent was  [-1,.5] the number of lattice points 
would be spaced according to 
\preformatted{
delta scalings: 0.6666667 0.3333333 0.1666667
}
Thus  \code{seq( -1,.5, 0.6666667 )} are the  (three) generated
points  within the spatial domain,  and of course 5 would also added beyond each endpoint.

To query the LKinfo object this information is in the 
\code{latticeInfo} component. E.g. the first set of lattice locations. 
\preformatted{
LKinfo$latticeInfo$grid[[1]]
$x
 [1] -4.3333 -3.6667 -3.0000 -2.3333 -1.6667
 [6] -1.0000 -0.3333  0.3333  1.0000  1.6667
[11]  2.3333  3.0000  3.6667  4.3333

$y
 [1] -4.3333 -3.6667 -3.0000 -2.3333 -1.6667
 [6] -1.0000 -0.3333  0.3333  1.0000  1.6667
[11]  2.3333  3.0000  3.6667  4.3333

attr(,"class")
[1] "gridList"
}

By default equal spacing of the lattice is assumed in the x and y
directions. To change this use the optional \code{V} argument to scale
the coordinates. For example, to set 4 lattice points in both
dimension for the above example:
\preformatted{
sDomain2<- cbind( c(-1,1), c( -1, .5))
LKinfo<- LKrigSetup(sDomain2, nlevel=3, NC=4,
                 a.wght=4.1, nu=1.0,
                 V= diag(c(2, 1.5)) )

print(LKinfo$latticeInfo$mx)

     [,1] [,2]
[1,]   14   14
[2,]   17   17
[3,]   23   23
}

\strong{About the basis functions} With the lattice points defined the
default basis functions are radial Wendland functions centered at
each point. The scaling of the basis functions is determined by the
\code{overlap} and the \code{delta} values.
In R code, if the lattice point is x0  at level l then the basis
function evaluated at location x1 is 
\preformatted{
basisFunctionValue <- Wendland( rdist( x1,x0)/( delta[l]*overlap))
}

If V is included then  x1 is transformed as  \code{x1\%*\%solve(V)}
before evaluating in the basis function. 

\strong{About  \code{a.wght}} 

For this geometry the basic form of awght is as a list with as many components as levels. However, the LKrigSetup function will reshape a scalar or vector argument in this this format. \code{a.wght} can take the following forms:

\strong{Scalar value: } In this case the value is used at all lattice points and at all levels in the SAR. Four nearest neighbors are set to -1. 

\strong{List/vector of length \code{nlevel} } In this case separate values for the a.wght 
will be used for the central SAR  value at each level.

\strong{List of vectors}  With a list of length \code{nlevel} and each
component is a 9 element vector, the values in the vector correspond to
the central lattice point and 8 nearest neighbors with the indexing:
\preformatted{
      1 4 7
      2 5 8
      3 6 9
}
E.g. the 5 element is the center, 3  is the lower left hand
corner etc. 

\strong{Non-stationary models}
This is the nuclear option to handle a completely non-stationary correlation
structure! In this case one specifies the models by passing prediction
objects for one or more of \code{rho.object}, \code{alphaObject} or 
\code{a.wghtObject}. 
}
\author{
Doug Nychka
%%  ~~who you are~~
}


\seealso{
\code{ \link{LatticeKrig}}
\code{ \link{LKGeometry}}
\code{ \link{LKrig}}

}
\examples{
# the grid with only 2 extra boundary points
  sDomain<- cbind( c(-1,1), c( -1, 1))
  LKinfo<- LKrigSetup(sDomain, nlevel=3, NC=4, a.wght=4.1,
           NC.buffer=2, alpha=c(1,.5,.125) )
  LKgrid<- LKinfo$latticeInfo$grid
  plot(   make.surface.grid(LKgrid[[1]]),
                          pch=16, cex=1.5)
  points( make.surface.grid(LKgrid[[2]]), 
                          pch=15, cex=.8, col="red" )
  points( make.surface.grid(LKgrid[[3]]),
                          pch="+", col="green" )
  rect(sDomain[1,1],sDomain[1,2],
     sDomain[2,1],sDomain[2,2], lwd=3 )

# basis functions on a grid
# this function actually evaluates all of them on the grid.
  xg<- make.surface.grid(
        list(x=seq( -2,2,,80), y=seq( -2,2,,80)) )
  out<- LKrig.basis( xg, LKinfo)
# basis functions 20, 26, 100  and 200
  plot(   make.surface.grid( LKgrid[[1]] ) , 
                          pch=16, cex=.5)
  rect(sDomain[1,1],sDomain[1,2],
     sDomain[2,1],sDomain[2,2], lwd=3,border="grey" )
  contour( as.surface(xg, out[,20]), col="red1",
                                     add=TRUE)
  contour( as.surface(xg, out[,36]), col="red4", 
                                     add=TRUE)
  contour( as.surface(xg, out[,100]), col="blue1",
                                     add=TRUE)
  contour( as.surface(xg, out[,200]), col="blue4",
                                      add=TRUE)
  title( "basis functions 20, 26, 100, 200")
}
\keyword{spatial}% __ONLY ONE__ keyword per line
