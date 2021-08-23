# LatticeKrig
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

suppressMessages(library(LatticeKrig))
 options( echo=FALSE)
 test.for.zero.flag<- 1
# 
# tests of LK computations with generalized least sqaures for the 
# fixed linear part of the model -- they should be identical!
#
  set.seed(123)
data(ozone2)

 x<-ozone2$lon.lat
 y<- ozone2$y[16,]
good<- !is.na( y)
x<- x[good,]
y<- y[good]

LKInfo<- LKrigSetup( x,y, NC=4, nlevel=2, a.wght=5, nu=1, lambda=.1 )
obj<- LKrig( x,y, LKinfo = LKInfo )
             
# create fixed part estimates using GLS formula
N<- length( y)

Sigma<- LKrig.cov( x,x, LKinfo=obj$LKinfo) +
  LKInfo$lambda  * diag( 1, length(y))

# classisc GLS by  inverse square root trick
X<- cbind( rep(1,length(y)), x)

temp<- eigen( Sigma, symmetric=TRUE)
U<- temp$vectors
D<- temp$values
# inverse square root of Sigma
Sinv2<-  U%*% diag( 1/sqrt(D))%*%t(U)
# transform  linear model to iid errors using inverse square root
#  proportional to covariance matrix
yStar<- Sinv2%*%y
XStar<- Sinv2%*%X

objlm<- lm( yStar~XStar - 1)
test.for.zero( coef(objlm), obj$d.coef)

# compare standard errors

covM0<- summary(objlm)$cov
SE0<- summary(objlm)$coefficients[,2]
lmSigma<- summary(objlm)$sigma
SE1<- sqrt(diag( covM0)) * summary(objlm)$sigma
covM1<- solve( t(XStar)%*%XStar) 
covM2<- obj$Omega
test.for.zero( SE0, SE1, tag="trans data to uncorrelated and GLS")
test.for.zero( covM0, covM1, 
               tag="cov matrix of parameters GLS")
test.for.zero( covM0, covM2,
               tag="cov matrix of parameters GLS from LatticeKrig")

SE2<- sqrt( diag( obj$rho.MLE* covM2) )

# checking how GLS residual SS is computed
SS0<- sum( objlm$residuals^2)

test.for.zero( obj$quad.form, SS0, tag="testing quadratic form")

test.for.zero( obj$rho.MLE*N, SS0, tag="testing rho*N and GLS SS")

test.for.zero( sqrt(obj$rho.MLE*(N/(N-3))),lmSigma , 
               tag="Estimates of rho and from GLS")

# now check with Z option do this by adding a quadratic part in Z
 allTerms<- fields.mkpoly(x, m=3)
 Z<- allTerms[,4:6]
 
 LKInfo<- LKrigSetup( x,y, NC=4, nlevel=2, a.wght=5, nu=1, 
                      lambda=.1, fixedFunctionArgs = list(m = 3) )
 objQuad<- LKrig( x,y, LKinfo = LKInfo )
 
 LKInfoZ<- LKrigSetup( x,y, NC=4, nlevel=2, a.wght=5, nu=1, 
                      lambda=.1, fixedFunctionArgs = list(m = 2) )
 objZ<- LKrig( x,y,Z=Z, LKinfo = LKInfoZ )
 test.for.zero( summary(objQuad)$coefficients, summary(objZ)$coefficients )
 
 X1<- cbind( X,Z)
 yStar<- Sinv2%*%y
 XStar<- Sinv2%*%X1
 
 objlm<- lm( yStar~XStar - 1)
 test.for.zero( coef(objlm), objQuad$d.coef)
 
 
 
 
 
 




