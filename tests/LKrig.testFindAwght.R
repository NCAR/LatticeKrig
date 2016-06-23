
# simulation example for estimating covariance
# LatticeKrig
# Copyright 2004-2016, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library(LatticeKrig)
options(echo = FALSE)
test.for.zero.flag <- 1

set.seed( 122)
M<- 1e2 # nummber of independent spatial replications
N<- 5e3 # number of obs
x<- matrix( runif(2*N ), N,2)
lambdaTrue<- .1^2
#NOTE: true rho is 1.0
LKinfoTrue<- LKrigSetup(x,NC=8, nlevel=3, a.wght= 4.2,
                        alpha=c(1.0,.5,.25), lambda=lambdaTrue,
                        normalize=FALSE, NC.buffer=2)
f<- LKrig.sim( x, LKinfoTrue, M=M) 
E<- matrix( rnorm( prod( dim( f))), nrow= nrow( f), ncol=ncol(f) )
Y<- f + sqrt(lambdaTrue)* E

LKinfoTest<- LKinfoUpdate( LKinfoTrue, a.wght=4.3,
                         lambda=lambdaTrue*(1.05))



Fit1<- LKrigFindLambdaAwght( x,Y,LKinfo=LKinfoTest,
                             verbose=FALSE)
print(Fit1$summary)

LKinfoTest2<- LKinfoUpdate(LKinfoTest, a.wght= Fit1$a.wght.MLE)
Fit2<- LKrigFindLambda( x,Y,LKinfo=LKinfoTest2)
Fit2$summary
Fit2$lambda.MLE

test.for.zero( Fit2$summary["lnProfLike"],
               Fit1$summary["lnProfLike"], tol=2e-7, tag=" lambda MLE")





