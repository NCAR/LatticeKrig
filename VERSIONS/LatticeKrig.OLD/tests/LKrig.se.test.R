# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library(LatticeKrig)
options( echo=FALSE)
 test.for.zero.flag<- 1

# small test dataset
data(ozone2)
x<- ozone2$lon.lat
y<- ozone2$y[16,]
good<- !is.na( y)
x<- x[good,]
y<- y[good]
x<- x[1:20,]
y<- y[1:20]

  a.wght<- 5
  lambda <-  1.5
# in both calls iseed is set the same so we can compare 
# Monte Carlo estimates of effective degrees of freedom
  obj1<- LKrig( x,y,NC=16,nlevel=1, alpha=1,  lambda=lambda, a.wght=5, NtrA=20,iseed=122,
               return.cholesky=TRUE)
  obj2<- mKrig( x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
                      cov.args=list( LKinfo=obj1$LKinfo), NtrA=20, iseed=122)
  obj3<- Krig( x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
                      cov.args=list( LKinfo=obj1$LKinfo) )
  test.for.zero( obj1$fitted.values, obj2$fitted.values,
                  tag="comparing predicted values LKrig and mKrig")
  test.for.zero( obj1$fitted.values, obj3$fitted.values,
                  tag="comparing predicted values LKrig and Krig")
             
  test3.se<-predictSE( obj3, x[1:3,])
  test2.se<-predictSE.mKrig( obj2, x[1:3,])
  test1.se<- predictSE.LKrig( obj1,  x[1:3,])
  test.for.zero( test1.se, test2.se,
                  tag="comparing SE values LKrig and mKrig equal weights")
  test.for.zero( test1.se, test3.se, 
                  tag="comparing SE values LKrig and Krig equal weights")

###########################
##### weighted case
###########################
 set.seed(123)
 weights<- runif( 20)
 obj1<- LKrig( x,y,NC=16,nlevel=1, alpha=1,  lambda=lambda, a.wght=5, NtrA=20,iseed=122,
               return.cholesky=TRUE, weights=weights) 
 obj2<- mKrig( x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
                      cov.args=list( LKinfo=obj1$LKinfo), NtrA=20, iseed=122,
                       weights=weights)
 obj0<- Krig(  x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
                      cov.args=list( LKinfo=obj1$LKinfo), weights=weights)
 test0.se<- predictSE.Krig( obj0, x[1:3,])
 test1.se<- predictSE.LKrig( obj1,  x[1:3,])
 test2.se<-predictSE.mKrig( obj2, x[1:3,])
 test.for.zero( test0.se,test1.se,
                  tag="sanity for SE Krig and mKrig unequal weights")
 test.for.zero( test1.se,test2.se,
                  tag=" SE LKrig and mKrig unequal weights")

###########################
######## covariates
##########################
 set.seed(122)
Z<- runif(20)
 obj1<- LKrig( x,y,NC=16,nlevel=1, alpha=1,  lambda=lambda, a.wght=5, NtrA=20,iseed=122,
               return.cholesky=TRUE, weights=weights, Z=Z)
 obj0<- Krig(  x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
                      cov.args=list( LKinfo=obj1$LKinfo), weights=weights, Z=Z)
  test0<- predictSE( obj0, drop.Z=FALSE, Z=Z)
  test1<- predictSE( obj1, drop.Z=FALSE, Z=Z)
  test.for.zero( test0, test1, tag="check on SE values with drop.Z=FALSE")
#
  test0<- predictSE( obj0, drop.Z=TRUE)
  test1<- predictSE( obj1, drop.Z=TRUE)
  test.for.zero( test0, test1, tag="check on SE values with drop.Z=TRUE")

cat("All done with SE tests", fil=TRUE)
options( echo=TRUE)

