
# simulation example for estimating covariance
# LatticeKrig
# Copyright 2004-2016, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library(LatticeKrig)
options(echo = FALSE)
test.for.zero.flag <- 1

set.seed( 122)


M<- 1e3 # nummber of independent spatial replications
N<- 1e2 # number of obs
x<- matrix( runif(2*N ), N,2)
lambdaTrue<- .1^2
#NOTE: true rho is 1.0 dont add fixed function so likelihood is precise.
LKinfoTrue<- LKrigSetup(x,NC=3, nlevel=3, a.wght= 4.2,
                        alpha=c(1.0,.5,.25), lambda=lambdaTrue,
                        normalize=FALSE, NC.buffer=2, fixedFunction = NULL)
f<- LKrig.sim( x, LKinfoTrue, M=M) 
E<- matrix( rnorm( prod( dim( f))), nrow= nrow( f), ncol=ncol(f) )
Y<- f + sqrt(lambdaTrue)* E

LKinfoTest<- LKinfoUpdate( LKinfoTrue, a.wght=4.3,
                         lambda=lambdaTrue*(1.05))


tick<- Sys.time()
Fit1<- LKrigFindLambdaAwght( x,Y,LKinfo=LKinfoTest,
                             verbose=FALSE)
tock<- Sys.time()
print( tock - tick)

print(Fit1$summary)

##plot( log10(Fit1$lnLike.eval[,1]), Fit1$lnLike.eval[,5], xlim=c(-3,-1.5), ylim=c(328000,339000))
#points( log10(Fit1$lambda.MLE),Fit1$summary["lnProfLike"], col="red")

#plot( log10(Fit2$lnLike.eval[,1]), Fit2$lnLike.eval[,4], xlim=c(-3,-1.5), ylim=c(328000,339000))
#points( log10(Fit2$lambda.MLE),Fit2$summary["lnProfLike"], col="red")

LKinfoTest2<- LKinfoUpdate(LKinfoTest, a.wght= Fit1$a.wght.MLE, lambda=Fit1$lambda.MLE )
Fit2<- LKrigFindLambda( x,Y,LKinfo=LKinfoTest2)
Fit2$summary
Fit2$lambda.MLE

test.for.zero( Fit2$summary["lnProfLike"],
               Fit1$summary["lnProfLike"], tol=1e-6, tag=" lambda MLE")





