# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html


# test of radial basis function based on Wendland
# and using sparse formats
# Important check is of the FORTRAN function dfind2d
# that does pairwise distances among points within a specified range.

  library(LatticeKrig)
  options( echo=FALSE)
  test.for.zero.flag<-1

  nLocation<- 4
  set.seed( 333)
  xLocation<- cbind( runif( 4, 3,5), runif( 4,3,5))
  LKinfo<- LKrig.setup( cbind( c(3,5), c(3,5)), NC=3, NC.buffer=2,
                                        a.wght=5,  alpha=1, nlevel=1, normalize=FALSE)
  testVar1<- LKrig.normalize.basis.fast(1, LKinfo, xLocation )
  testVar0<-LKrig.cov( xLocation,   LKinfo= LKinfo, marginal=TRUE)
  test.for.zero( testVar0, testVar1, tag="Marginal variance and fast normalize")
###  multiple levels

  alphaVec<- c( 1,.5,.2)
  LKinfo<- LKrig.setup( cbind( c(3,5), c(3,5)), NC=3, NC.buffer=2,
                                        a.wght=5,  alpha=alphaVec, nlevel=3, normalize=FALSE)
  test1<- LKrig.normalize.basis.fast(1, LKinfo, xLocation )
  test2<- LKrig.normalize.basis.fast(2, LKinfo, xLocation )
  test3<- LKrig.normalize.basis.fast(3, LKinfo, xLocation )
  testVar1<- cbind( test1, test2, test3) %*% alphaVec
  testVar0<-LKrig.cov( xLocation,   LKinfo= LKinfo, marginal=TRUE)
  test.for.zero( testVar0, testVar1, tag="Marginal variance and fast normalize")

cat( "Done with testing fast normalize algorithm", fill=TRUE)
options( echo=TRUE)

