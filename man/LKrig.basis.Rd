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

\name{LKrig.basis}
\alias{LKrig.basis}
\alias{LKrig.precision}
\alias{LKrig.cov}
\alias{LKrig.cov.plot}
\alias{LKrig.quadraticform}
\alias{LKrig.spind2spam}
\alias{LKrigCovWeightedObs}
\alias{LKrigMarginalVariance}
\title{
Functions for generating a multiresolution, compactly supported basis,  
multiresolution covariance functions and simulating from these processes.
}
\description{
These functions support the \code{LKrig} function.  Their main
function is to create and evaluate radial basis functions of varying
support on a nested set of regular grids.  This series of grids forms
a multiresolution basis.  The Gaussian process model is an expansion
in these basis functions where the basis coefficients follow a Markov
random field model for each resolution level.  This family of
functions generate the basis using sparse matrices, evaluate the
covariance function of the process and also simulate realizations of
the process. \code{LKrig.cov.plot} is a useful function to get a quick
plot of the covariance function implied by a LatticeKrig
specification.}

\usage{
#
LKrig.cov(x1, x2 = NULL, LKinfo, C = NA, marginal = FALSE) 
LKrig.cov.plot( LKinfo, NP=200, center = NULL, xlim = NULL, ylim = NULL)
LKrigCovWeightedObs(x1, wX, LKinfo) 
LKrig.basis(x1, LKinfo, verbose = FALSE)
LKrig.precision(LKinfo, return.B = FALSE, verbose=FALSE)  
LKrig.quadraticform( Q, PHI, choleskyMemory = NULL) 
LKrig.spind2spam(obj, add.zero.rows=TRUE)
LKrigMarginalVariance(x1, LKinfo, verbose = FALSE)
}

\arguments{

\item{add.zero.rows}{If TRUE the conversion of the sparse matrix to
spam format will have at least one element in each row. If there are
no elements explicitly given in \code{obj} then an element with value
zero is added. This technical detail is needed to accommodate the spam
format for sparse matrices.}


\item{C}{If passed the covariance matrix will be multiplied by this
vector or matrix.}

\item{center}{The point in the spatial domain that is used to evaluate the
covariance function. The evaluation is done on x and y transects
through the spatial domain intersecting at \code{center} and finding
the covariance with respect to this point. If NULL defaults to the
center of the spatial domain.}

\item{choleskyMemory}{A list giving the memory requirements for 
	a sparse cholesky decomposition. See \code{chol.spam} for details and
	also \code{LKrig}. }
	

\item{LKinfo}{A list with components that give the information
describing a multiresolution basis with a Markov random field used for
the covariance of the basis coefficients. This list is created in
\code{LKrig} or by \code{LKrigSetup} and returned in the output
object.  (See section on returned Value below for this list's
description.)}

\item{marginal}{If TRUE returns the marginal variance. Currently not
implemented!}

\item{NP}{Number of points to evaluate the covariance function along
each transect of the spatial domain.}

\item{obj}{An object returned by \code{LKrig} or a sparse matrix in
row/column format passed to \code{LKrig.spind2spam}.}

\item{PHI}{ A sparse matrix of basis functions (rows index points for
evaluation, columns index basis functions).}

\item{return.B}{ If TRUE B is returned instead of the precision matrix
t(B)\%*\%B.}

\item{Q}{A sparse (spam format) precision matrix.}

\item{wX}{The wX matrix constructed in the \code{\link{LKrig}} function.}

\item{x1}{A two column matrix of 2-dimension locations to evaluate
basis functions or the first set of locations to evaluate the
covariance function or the locations for the simulated process.  Rows
index the different locations: to be precise \code{x1[i,1:2]} are the
"x" and "y" coordinates for the \code{i} th location. }

\item{x2}{Second set of locations to evaluate covariance function.}

\item{xlim}{Limits in x coordinate for evaluating the covariance
model. Default is the spatial domain.}

\item{ylim}{Limits in y coordinate for evaluating the covariance
model.Default is the spatial domain.}

\item{verbose}{If TRUE intermediate steps and other debugging
information are printed.}

} 
\details{ The basis functions are two-dimensional radial basis
functions based on the compactly supported stationary covariance function
(Wendland covariance) and centered on
regular grid points with the scaling tied to the grid spacing.

For a basis at the coarsest level, the grid centers are generated by
expanding the two sequences
\preformatted{
seq(grid.info$xmin,grid.info$xmax,grid.info$delta)
seq(grid.info$ymin,grid.info$ymax,grid.info$delta) 
}
into a regular grid of center points. The same spacing \code{delta} is
used in both directions.  The unnormalized basis functions are
evaluated at locations \code{x1} by finding the pairwise, radial
distances among \code{centers} and \code{x1}, scaling by
\code{grid.info$delta * overlap} and then evaluating with the function
name passed as \code{BasisFunction}. By default this is the 2-d
Wendland covariance of order 2.  Perhaps the most important point
about the \code{LKrig.basis} is that it is designed to return a matrix
of all basis functions as a sequence of points. There is no need to
have a function that evaluates basis functions individually.  In R
code for a set of locations \code{x1} and a rectangular spatial domain
with ranges \code{xmin, xmax, ymin ,ymax}:

\preformatted{
 centers<- expand.grid(seq(xmin,xmax,delta),
                       seq(ymin,ymax,delta) )
 bigD<- rdist( x1, centers)/(delta*2.5)
 PHI<- Wendland.function( bigD)
} 

Note that there will be \code{nrow(centers)} basis functions generated
where the precise number depends on the range of the domain and the
choice of delta. The basis functions are indexed by the columns in
\code{PHI} and this is a convention throughout this package.  There
will be \code{nrow(x1)} rows in \code{PHI} as each basis function will be
evaluated at each 2-d location.

The basis functions are then normalized by scaling the basis functions
at each location so that resulting marginal variance of the process is
1. This is done to coax the covariance model closer to a stationary
representation. It is easiest to express this normalization by pseudo
R code:

If \code{Q} is the precision matrix of the basis coefficients then 
in R/LatticeKrig code:
\preformatted{
Omega<-  solve(Q)
process.variance <- diag(PHI\%*\% Omega \%*\%t(PHI) )
PHI.normalized <-  diag(1/sqrt(process.variance)) \%*\% PHI
}
where \code{Omega} is the unnormalized covariance matrix of the basis
function coefficients.

Although correct, the code above is not an efficient algorithm to
compute the unnormalized process variance.  First the normalization
can be done level by level rather than dealing with the entire
multiresolution process at once. Also it is important to work with the
precision matrix rather than the covariance.  The function
\code{LKrig.normalize.basis.slow} takes advantage of the sparsity of the
precision matrix for the coefficients.
\code{LKrig.normalize.basis.fast} is even a more efficient version when
a.wght is constant for a given level and takes advantage of the
Kronecker structure on the precision matrix at each level.

The precision matrix for the basis coefficients at each resolution has
 the form \code{t(B)\%*\% B}. These matrices for the individual levels
 are assembled by \code{LKrig.precision} as the block diagonals of a
 larger precision matrix for the entire vector of coefficients. Note
 these matrices are created in a sparse format.  The specific
 entries in B, the object created by \code{LKrig.MRF.precision}, are a
 first order Markov random field: without edge adjustments the
 diagonal elements have the value \code{a.wght} and the first order
 neighbors have the value -1.

Below we give more details on how the weights are determined.
Following the notation in Lindgren and Rue a.wght= 4 + k2
with k2 greater than or equal to 0. Some schematics for filling
in the B matrix are given below (values are weights for the SAR on the
lattice with a period indicating zero weights).

\preformatted{
                                                 __________ 
   .   -1     .         |  -1      .            |  a.wght  -1
                        |                       |
  -1  a.wght  -1        | a.wght   -1           |   -1      .
                        |                       |
   .  -1      .         |  -1      .            |    .      .

Interior point        Left edge              Upper left corner
}

Previous versions of LatticeKrig considered an edge correction to
reflect other boundary conditions.  We have found these corrections to
be numerically unstable, however, and so prefer at this time of
writing adding a buffer of lattice points and using the uncorrected
weights described above.
 
}

\value{

\strong{LKrig.basis:} A matrix with number of rows equal to the rows of
\code{x1} and columns equal to the number of basis functions
(LKinfo$m). Attached to the matrix is an \code{info} attribute that
contains the list passed as \code{LKinfo}. Usually this value is in
spam sparse matrix format.


\strong{LKrig.precision:} For \code{return.B ==FALSE} a sparse, square
matrix with dimensions of the number of basis functions. For
\code{return.B == TRUE} the "B" SAR matrix is returned. This is useful
for checking this function.


\strong{LKrig.cov:} If \code{C=NA} a cross covariance matrix
with dimensions \code{nrow(x1)} and \code{nrow(x2)} is used. If \code{C} is
passed the result of multiplying the cross covariance matrix times
\code{C} is used.

\strong{LKrigMarginalVariance:} Gives the marginal variance of the
LatticeKrig process at each level and at the locations in \code{x1}. Returned value is a matrix with
\code{nlevel} columns indexing the levels and the number of rows equal to  \code{nrow(x1)}. 
If  \code{varLevels} is a row of this matrix then  sum( varLevels* LKinfo$alpha) is the marginal variance
for the full process when the different levels are weighted by alpha. This is weighted sum is the marginal
variance returned by \code{LKrig.cov} and \code{ marginal=TRUE} ( Also assuming that \code{LKinfo$rho.object}
is NULL, which it usaully is.) . 




\strong{LKrig.sim:} A matrix with dimensions of \code{nrow(x1)} by
\code{M}. Each column are vectors of simulated values at the locations
\code{x1}.

\strong{LKrig.cov.plot:} Evaluates the covariance specified in the list
LKinfo with respect to the point \code{center} along a transects in
the x and y directions intersecting this point.  Note the rectangular
extent of the spatial domain is part of the grid information in
LKinfo.  Returns components \code{u}, \code{d} and \code{cov}. Each of
these are two column matrices with the first column being results in
the x direction and second column in the y direction. \code{d} are the
distances of the points from the center and \code{u} are the actual x
or y spatial coordinates. \code{cov} are the values of the covariance
function. If normalize is TRUE these will in fact be the correlation
functions.  To plot the returned list use
\preformatted{
 out<- LKrig.cov.plot(LKinfo)
 matplot( out$d, out$cov, type="l")
}

\strong{LKrig.quadraticform:} Returns a vector that is
 \code{diag(t(PHI)\%*\% solve(Q) \%*\% PHI))} closely related to the marginal
 variance of the process.

\strong{LKrig.normalize.basis},
\strong{LKrig.normalize.basis.fast}, 
\strong{LKrig.normalize.basis.slow}: A vector of variances
corresponding to the unnormalized process at the locations.


\strong{LKrig.spind2spam:} This converts a matrix in spind sparse format to spam sparse format. Although useful for checking, LatticeKrig now uses the
the \code{spam} function directly to do the conversion. 
If obj is a list in spind format then the 

\code{ obj2<- LKrig.spind2spam(obj)} 

and

\code{ obj2<- spam(obj[c("ind","ra")], nrow=obj$da[1], ncol=obj$da[2])}

give the same result. 

}

\author{
Doug Nychka
}

\seealso{
LKrig, mKrig, Krig, fastTps, Wendland, LKrigSAR
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

# What does the  LatticeKrig covariance function look like?
# set up LKinfo object
# NC=10 sets the grid for the first level of basis functions
# NC^2 = 100 grid points in first level if square domain.
# given four levels the number of basis functions
# = 10^2 + 19^2 +37^2 + 73^2 = 5329
# effective range scales as roughly kappa where a.wght =  4 + kappa^2    
# or exponential decreasing marginal variances for the components.
    NC<- 10 
    nlevel<- 4
    a.wght<- rep(  4 + 1/(.5)^2, nlevel)
    alpha<-  1/2^(0:(nlevel-1)) 
    LKinfo2<- LKrigSetup( cbind( c( -1,1), c(-1,1)), NC=NC,
                   nlevel=nlevel, a.wght=a.wght,alpha=alpha)
# evaluate covariance  along the  horizontal line through
# midpoint of region -- (0,0) in this case. 
    look<- LKrig.cov.plot( LKinfo2)
# a plot of the covariance function in x and y with respect to (0,0)
    set.panel(2,1)  
    plot(look$u[,1], look$cov[,1], type="l")
    title("X transect")
    plot(look$u[,2], look$cov[,2], type="l")
    title("Y transect")
    set.panel(1,1)
#
#
\dontrun{
# full 2-d view of the covariance (this example follows the code
# in LKrig.cov.plot)
 x2<- cbind( 0,0)
 x1<- make.surface.grid( list(x=seq( -1,1,,40),  y=seq( -1,1,,40)))
 look<- LKrig.cov( x1,x2, LKinfo2)
 contour( as.surface( x1, look))
# Note nearly circular contours.
# of  course  plot(look[,80/2]) should look like plot above.
#
}

\dontrun{
#Some correlation functions from different models
set.panel(2,1)
# a selection of ranges:
  hold<- matrix( NA, nrow=150, ncol=4)
  kappa<- seq( .25,1,,4)
  x2<- cbind( 0,0)
  x1<-  cbind( seq(-1,1,,150), rep( 0,150))
  for( k in 1:4){
    LKtemp<-  LKrigSetup( cbind( c( -1,1), c(-1,1)), NC=NC,
                   nlevel=nlevel,
                   a.wght= 4  + 1/(kappa[k]^2),
                   alpha=alpha)
    hold[,k]<-  LKrig.cov( x1,x2, LKinfo=LKtemp)
  }
  matplot( x1[,1], hold, type="l", lty=1, col=rainbow(5), pch=16 )
# a selection of smoothness parameters
  ktemp<- .5 # fix range
  alpha.power<- seq( 1,4,,4)
  LKtemp<- LKinfo2
  for( k in 1:4){
   LKtemp<-  LKrigSetup( cbind( c( -1,1), c(-1,1)), NC=NC,
                   nlevel=nlevel,
                   a.wght= 4  + 1/(ktemp^2),
                   alpha=alpha^alpha.power[k])
    hold[,k]<-  LKrig.cov( x1,x2, LKinfo=LKtemp)
  }
  matplot( x1[,1], hold, type="l", lty=1, col=rainbow(5) )
 set.panel()
}
 
\dontrun{
# generating a basis on the domain [-1,1] by [-1,1] with 1 level
# Default number of buffer points are added to each side. 
  LKinfo<- LKrigSetup(cbind( c(-1,1), c(-1,1)), NC=6,
                                 nlevel=1, a.wght=4.5,alpha=1, NC.buffer=0 )
# evaluate the basis functions on a grid to look at them
  xg<- make.surface.grid( list(x=seq(-1,1,,50), y= seq(-1,1,,50)))
  PHI<- LKrig.basis( xg,LKinfo)
  dim(PHI) # should be  2500=50^2  by  36=6^2
# plot the 9th basis function  as.surface is a handy function to
# reformat the vector as an image object
# using the grid information in an attribute of the grid points
  set.panel(1,3) 
  image.plot(as.surface(xg, PHI[,9]))
  points(  make.surface.grid( LKrigLatticeCenters(LKinfo, 1)) , col="grey", cex=.5)
  title("A radial basis function")
# compare to the tensor product basis type
  LKinfo2<- LKrigSetup(cbind( c(-1,1), c(-1,1)), NC=6,
                                 nlevel=1, a.wght=4.5,alpha=1, NC.buffer=0,
                                 BasisType="Tensor" )
  PHI2<- LKrig.basis( xg,LKinfo2)
  image.plot(as.surface(xg, PHI2[,9]))
  points(  make.surface.grid( LKrigLatticeCenters(LKinfo, 1)), col="grey", cex=.5)
  title("Tensor product basis function")
  
  image.plot(as.surface(xg, PHI[,9] - PHI2[,9]))
  points(  make.surface.grid( LKrigLatticeCenters(LKinfo, 1)), col="grey", cex=.5)
  title(" Radial - Tensor for 9th basis function")                       
set.panel()
}
#
# example of basis function indexing
#
\dontrun{
# generating a basis on the domain [-1,1]X[-1,1] with 3 levels
# note that there are no buffering grid points.
  set.panel(3,2)
  LKinfo<-LKrigSetup(cbind( c(-1,1), c(-1,1)), NC=6,
                    a.wght=rep(5,3), alpha=c(1,.5,.25), nlevel=3,
                    NC.buffer=0)
# evaluate the basis functions on a grid to look at them
  xtemp<- seq(-1,1,,40)
  xg<- make.surface.grid( list(x=xtemp, y= xtemp) )
  PHI<- LKrig.basis( xg,LKinfo)
# coerce to dense matrix format to make plotting easier.
  PHI<- spam2full(PHI)
# first tenth, and last basis function in each resolution level
# basis functions centers are added
 set.panel(3,3)
    for(  j in 1:3){
      id1<- LKinfo$latticeInfo$offset[j]+ 1
      id2<-  LKinfo$latticeInfo$offset[j]+ 10
      idlast<- LKinfo$latticeInfo$offset[j] +
                  LKinfo$latticeInfo$mx[j,1]*LKinfo$latticeInfo$mx[j,2]
   
      centers<-  make.surface.grid(LKrigLatticeCenters(LKinfo, j) )
      image.plot( as.surface(xg, PHI[,id1]))
      points( centers, cex=.2, col="grey")
      image.plot(as.surface(xg, PHI[,id2]))
      points( centers, cex=.2, col="grey")
      image.plot( as.surface(xg, PHI[,idlast]))
      points( centers, cex=.2, col="grey")
}
  set.panel()
}
\dontrun{
# examining the stationarity of covariance model
  temp.fun<- 
     function( NC.buffer=0, NC=4,  a.wght=4.01){
        LKinfo<- LKrigSetup(cbind( c(-1,1), c(-1,1)),nlevel=1, alpha=1,
                                 a.wght=a.wght, NC=NC,   
                                 NC.buffer=NC.buffer,
                                  choleskyMemory=list(nnzR=2e6))
        cov1y<- cov1x<- cov0x<- cov0y<-  matrix( NA, nrow=200, ncol=20)
        cov1dx<- cov1dy<- cov0dx<- cov0dy<- matrix( NA, nrow=200, ncol=20)
        cgrid<- seq( 0,1,,20)
        for( k in 1:20){
            hold<- LKrig.cov.plot( LKinfo,
                            center=rbind( c(cgrid[k], cgrid[k])), NP=200)
            cov1x[,k] <- hold$cov[,1]
            cov1y[,k] <- hold$cov[,2]
            cov1dx[,k] <- hold$d[,1]
            cov1dy[,k] <- hold$d[,2]
            hold<- LKrig.cov.plot( LKinfo,
                             center=rbind( c(cgrid[k],0) ), NP=200)
            cov0x[,k] <- hold$cov[,1]
            cov0y[,k] <- hold$cov[,2]
            cov0dx[,k] <- hold$d[,1]
            cov0dy[,k] <- hold$d[,2]
                }
         matplot( cov1dx, cov1x, type="l", col= rainbow(20),
                         xlab="", ylab="correlation")
         mtext( side=1, line=-1, text="diagonal X")
         title( paste(  " buffer=",NC.buffer), cex=.5)
         matplot( cov1dy, cov1y, type="l", col= rainbow(20),
                        xlab="", ylab="")
         mtext( side=1, line=-1, text="diagonal Y")
         matplot(cov0dx, cov0x, type="l", col= rainbow(20),
                        xlab="",       ylab="")
         mtext( side=1, line=-1, text="middle X")
         matplot( cov0dy, cov0y, type="l", col= rainbow(20),
                         xlab="",   ylab="")
         mtext( side=1, line=-1, text="middle Y")
         title( paste( NC, a.wght), cex=.5)
}


 set.panel(3,4)
par(mar=c(3,4,1,0), oma=c(1,1,1,1))
temp.fun(  NC.buffer=5, NC=4, a.wght=4.05)
temp.fun(  NC.buffer=5, NC=16, a.wght=4.05)
temp.fun(  NC.buffer=5, NC=64, a.wght=4.05)

set.panel(4,4)
par(mar=c(3,4,1,0), oma=c(1,1,1,1))
temp.fun( NC.buffer=0, NC=8)
temp.fun( NC.buffer=2, NC=8)
temp.fun( NC.buffer=4, NC=8)
# this next one takes a while
temp.fun( NC.buffer=8,  NC=8)
# stationary == curves should all line up!

}

}

\keyword{spatial}


