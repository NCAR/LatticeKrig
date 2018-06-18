# LatticeKrig
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html


# test of radial basis function based on Wendland
# and using sparse formats
# Important check is of the FORTRAN function dfind2d
# that does pairwise distances among points within a specified range.

library(LatticeKrig)
options( echo=FALSE)
 test.for.zero.flag<- 1
set.seed(123)
x1<-  matrix( runif(40), ncol=2)
x2<-  matrix( runif(30), ncol=2)
n1<- nrow(x1)
n2<- nrow(x2)

look1<-Radial.basis(x1,x2, delta=.7, just.distance=TRUE)
look1<- spind2full(look1)
look2<- rdist( x1,x2)
look2[ look2>.7] <-0
test.for.zero( look1,look2)

# test when range varies among different points
delta<- c( rep(.6,10), rep( .3,n2-10))
look1<- Radial.basis(x1,x2, delta=delta,just.distance=TRUE)
look1<- spind2full(look1)
ind<-matrix( delta, nrow=n1,ncol=n2, byrow=TRUE)
look2<- rdist( x1,x2)
look2[ look2> ind] <- 0
test.for.zero( look1,look2)

look1<-Radial.basis(x1,x2, delta=.5)
look1<- spam2full(look1)
look2<- Wendland2.2(rdist( x1,x2)/.5)
test.for.zero( look1,look2)

delta<- c( rep(.6,10), rep( .3,n2-10))
look1<-Radial.basis(x1,x2, delta=delta)
look1<- spam2full(look1)
dind<-matrix( delta, nrow=n1,ncol=n2, byrow=TRUE)
look2<- Wendland2.2(rdist( x1,x2)/dind)
test.for.zero( look1,look2)

# tests for cylindrical geometry
  x1<- make.surface.grid( list( x= seq(0,360,,120), y= seq( -30,30,,20) ))
  center<- rbind(c(20,5))
  R<- 360/(2*pi)
  temp<- rdist( LKrig.cyl(x1,R), (LKrig.cyl(center, R)))
  temp[ temp>=60] <- 0
#  image.plot( as.surface(x1,c(temp)))
  look1<-Radial.basis(x1,center, delta=60, just.distance=TRUE, distance.type="cylinder")
  look2<- spind2full(look1)
  test.for.zero( c(look2), c(temp), tag="cyl distance1")

  x1<- make.surface.grid( list( x= seq(0,360,,100), y= seq( -30,30,,15) ))
  center<-  make.surface.grid( list( x= seq(0,360,,15), y= seq( -30,30,,5) ))
  look1<-Radial.basis(x1,center, delta=60, distance.type="cylinder")
  look2<- spam2full(look1)
  R<- 360/(2*pi)
  temp<- rdist( LKrig.cyl(x1,R), (LKrig.cyl(center, R)))
  temp2<- WendlandFunction( temp/60)
  test.for.zero( c(look2), c(temp2), tag="cyl distance2")

  obj<- LKrig.setup( rbind( c(0,-30), c(360,30)), NC=10, a.wght=rep(5,4), nlevel=4,
                                  alpha=rep(1,4),distance.type="cylinder")


cat( "Done with testing LKrig basis", fill=TRUE)
options( echo=TRUE)
