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
% Earlier modifications DH 2/03/2012

\name{LKrig.MLE}
\alias{LKrig.MLE}
\alias{LKrigFindLambda}
\alias{LKrigFindLambdaAwght}
\alias{LambdaAwghtObjectiveFunction}
\alias{LKrig.make.par.grid}
\alias{omega2Awght}
\alias{Awght2Omega}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simple function to search over covariance parameters for Lattice Krig.}
\description{
Given a list of different covariance parameters for the Lattice Krig
covariance model this function computes the likelihood or a profiled
version (over lambda) and approximates a generalized cross-validation
function at each of the parameter settings. This is an experimental
function that has been productively used with a Latin hypercube design
package to efficiently search through the LatticeKrig covariance
parameter space.
}
\usage{
LKrigFindLambda(x, y, ..., LKinfo,
                            use.cholesky=NULL,                          
                            lambda.profile=TRUE,
                            lowerBoundLogLambda=-16,tol = 0.005,
                            verbose=FALSE)
			    
LKrigFindLambdaAwght(x, y, ..., LKinfo, use.cholesky = NULL,
                 lowerBoundLogLambda = -16, upperBoundLogLambda = 4,
                 lowerBoundOmega = -3, upperBoundOmega = 0.75, factr =
                 1e+07, pgtol=1e-1, maxit = 15, verbose = FALSE)

                                 
LambdaAwghtObjectiveFunction(PARS, LKrigArgs, capture.env, verbose=FALSE )

LKrig.MLE( x,y,..., LKinfo, use.cholesky = NULL,
                    par.grid=NULL,
                    lambda.profile=TRUE,
                    verbose=FALSE,
                    lowerBoundLogLambda = -16,
                    nTasks = 1, taskID = 1,
                    tol = 0.005)
LKrig.make.par.grid(par.grid=NULL, LKinfo = NULL) 
omega2Awght (omega, LKinfo)
Awght2Omega (Awght, LKinfo)
}
\arguments{
\item{Awght}{Value of Awght parameter to convert to omega form}
\item{capture.env}{The enviornment to save to the likelihood evaluation to the 
object \code{capture.evaluations}. }


\item{lambda.profile}{ A logical value controlling whether the
 likelihood is maximized over lambda.  For LKrigFindLambda if TRUE the
 likelihood is maximized over lambda at the covariance values in
 LKinfo and if FALSE the likelihood is just evaluated at LKinfo
 including the lambda value in this list.  For LKrig.MLE if TRUE for
 each set of parameters in par.grid the value of lambda is found that
 maximizes the likelihood.  In this case the llambda value is the
 starting value for the optimizer. If llammbda[k] is NA then the
 lambda value found from the k-1 maximization is used as a starting
 value for the k step. (In the source code this is
 \code{llambda.opt}.) Of course this only makes sense if the other
 parameters are ordered so that the results for k-1 make sense as a
 lambda starting value for k.  If FALSE the likelihood is evaluated
 for the covariance parameters at the kth positions in the par.grid
 list including lambda.}

\item{LKinfo}{ An LKinfo object that specifies the LatticeKrig
covariance. Usually this is obtained by a call to \code{LKrigSetup}
or as the component LKinfo from the LKrig object.  The search
sequentially replaces the alpha and a.wght arguments in this list by
the values in par.grid but leaves everything else the same. If
par.grid is not passed the parameter values in LKinfo are used to
evaluate the likelihood. This option is most useful if one has fixed
values of alpha and a.wght and the goal is to maximize the likelihood
over lambda.  }
\item{LKrigArgs}{Argument list to call LKrig.}

\item{lowerBoundLogLambda}{Lower limit for lambda in searching for
MLE.}
\item{lowerBoundOmega}{Lower limit for omega in searching for
MLE.}
\item{maxit}{Maximum number of iterations (passed to optim function)}
\item{nTasks}{If using Rmpi the number of slaves available.}
\item{omega}{Value of the omega parameter to convert to  the a.wght
format.  }
\item{PARS}{ PAR[1]= log(lambda) and PAR[2]= .5*log(a.wght-4)
(also referred to as \code{omega}).
For the LKRectangle
geometry.}
 \item{par.grid}{ A list with components llambda, alpha, a.wght giving
the different sets of parameters to evaluate. If M is the number of
parameter setting to evaluate llambda is a vector length M and alpha
and a.wght are matrices with M rows and nlevel columns.  Thus, the kth
trial has parameters par.grid\$llambda[k], par.grid\$alpha[k,] and
par.grid\$a.wght[k,]. Currently this function does not support passing
a non-stationary spatial parameterization for alpha. The LKinfo object
details the other parts of the covariance specification (e.g. number
of levels, grid sizes) that do not change.  Note that par.grid assumes
\emph{ln} lambda not lambda. See details below for some other features
of the par.grid arguments.  }
\item{factr}{Controls convergence for the BFGS-L method.
( passed to optim).}
\item{pgtol}{Mysterious tolerance for gradient converenge in L-BFGS. This
seems to influence the number of iterations the most.}
 
\item{tol}{Tolerance on log likelihood use to determine convergence.}

\item{taskID}{If using Rmpi the slave id.}

\item{verbose}{If TRUE prints out intermediate results.}
\item{upperBoundOmega}{Upper limit for omega in searching for
MLE.}
\item{upperBoundLogLambda}{Upper limit for log lambda in searching for
MLE.}
\item{use.cholesky}{If not NULL then this object is used as the
symbolic cholesky decomposition of the covariance matrix for computing
the likelihood.}
\item{x}{ The spatial locations.  }
\item{y}{ The observations.  }
\item{\dots}{Any arguments to be passed to LKrig. E.g. \code{x},
       \code{y}, \code{Z} a covariate or \code{weights}.}  }
\details{

\strong{LKrigFindLambda:} Uses a simple one dimensional optimizer 
\code{optimize}. To maximize the log likelihood for log lambda over the
range: llambda.start + [-8,5]. This function is used to determine lambda in 
\code{LatticeKrig}.

\strong{LKrigFindLambdaAwght:} Uses a simple optimizer 
\code{optim}. To maximize the log likelihood for  lambda  and  a.wght over the
range. 

\strong{LKrig.MLE:} This is a simple wrapper function to accomplish
repeated calls to the LKrig function to evaluate the profile
likelihood and/or to optimize the likelihood over the lambda
parameters. The main point is that maximization over the lambda
parameter (or equivalently for sigma and rho) is the most important
and should be done before considering variation of other parameters.  If
lambda is specified then one has closed form expressions for sigma,
rho that can then be substituted back into the log full
likelihood. This operation that is the default throughout LatticeKrig
(and fields) can concentrate the likelihood on a reduced set
of parameters. The further refinement when \code{lambda.profile==TRUE}
is to maximize the concentrated likelihood over lambda and report
this result. This will be a profile likelihood over the remaining
parameters of the covariance.

The covariance/model parameters are alpha, a.wght, and log lambda and
are separate matrix or vector components of the \code{par.grid}
list. The cleanest version of this function would just require the
\code{par.grid} list, however, to be easier to use there are several
options to give partial information and let the function itself create
the master parameter list. For example, just a search over lambda
should be easy and not require creating \code{par.grid} outside the
function. To follow this option one can just give an LKinfo
object. The value for the lambda component in this object will be the
starting value with the default starting value being lambda =1.0.


In the second example below most of the coding is getting the grid of
parameters to search in the right form.  It is useful to normalize the
alpha parameters to sum to one so that the marginal variance of the
process is only parameterized by rho.  To make this easy to implement
there is the option to specify the alpha parameters in the form of a
mixture model so that the components are positive and add to one (the
gamma variable below). If a component \code{gamma} is passed as a
component of \code{par.grid} then this is assumed to be in the mixture
model form and the alpha weights are computed from this. Note that
\code{gamma} will be a matrix with (nlevel - 1) columns while
\code{alpha} has nlevel columns.


For those readers that use \code{which.max} these functions are natural
extensions and are handy for looking at interpolated surfaces of the
likelihood function.

\strong{which.max.matrix:} Finds the maximum value in a matrix and
returns the row/column index.

\strong{ \code{which.max.image }} Finds the maximum value in an image matrix and
returns the index and the corresponding grid values.

\strong{LKrig.make.par.grid:} This is usually used as an internal
function that converts the list of parameters in par.grid and the
LKinfo object into an more complex data structure used by
LKrig.MLE. Its returned value is a "list of lists" to make the search
over different parameters combinations simple. 

\strong{omega2Awght, Awght2omega } Converts between the Awght parameter (the diagonal elements
of the SAR matrix) and the omega parameter that provides an unconstrained range for optimization. 
The link is \code{Awght <- LKinfo$floorAwght + exp(omega) * (xDimension)}. Note that LKinfo
supplied the lower bound on the Awght becasue this is geometry/problem specific. For the 2-d rectangle
this is 4. 

}
\value{

\strong{LKrigFindLambda} 

\item{summary}{Giving information on the optimization over lambda.}

\item{LKinfo}{Covariance information object.}

\item{llambda.start, lambda.MLE}{Initial and final values for lambda.} 

\item{lnLike.eval}{Matrix with all values of log likelihood that were evaluated}

\item{call}{Calling arguments.}

\item{Mc}{Cholesky decomposition.}

\strong{LKrig.MLE}  

\item{summary}{A matrix with columns: effective degrees of freedom, ln
Profile likelihood, Generalized cross-validation function, MLE sigma,
MLE rho, full likelihood and number of parameter evaluations. The
rows correspond to the different parameters in the rows of the
\code{par.grid} components. }

\item{par.grid }{List of parameters used in search. Some parameters
might be filled in from the initial par.grid list passed and also from
\code{LKinfo}.}  

\item{LKinfo}{ LKinfo list that was either passed or
created.}

\item{index.MLE}{Index for case that has largest Likelihood value.}

\item{index.GCV}{Index for case that has largest GCV value.}

\item{LKinfo.MLE}{LKinfo list at the parameters with largest profile
likelihood.}

\item{lambda.MLE}{Value of lambda from grid with largest profile
likelihood.}

\item{call}{Calling sequence for this function.}

\strong{which.max.matrix} Returns a 2 column matrix with row and column
index of maximum. 

\strong{which.max.image} For an object in image format 
returns components x,y,z giving the location and the maximum value
for the image.
Also included is the component ind that is  the row and column indices for the maximum in the image matrix. 

\strong{LKrig.make.par.grid} Returns a list with components, alpha,
a.wght. Each component is a list where each component of the list is a
separate set of parameters. This more general format is useful for the
nonstationary case when the parameters alpha might be a list of nlevel
matrices.

}

\author{Douglas Nychka}

\seealso{
\code{\link{LKrig} } \code{\link{LatticeKrig} }
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# 
# fitting summer precip for  sub region of North America
# total precip in 1/10 mm for JJA 
  data(NorthAmericanRainfall)
# rename for less typing
  x<- cbind( NorthAmericanRainfall$longitude, NorthAmericanRainfall$latitude)
  y<- log10(NorthAmericanRainfall$precip)
# cut down the size of this data set so examples run quickly
  ind<- x[,1] > -90 & x[,2] < 35 #
  x<- x[ind,]
  y<- y[ind]

# This is a single level smoother
 
  LKinfo<- LKrigSetup(x,NC=25,nlevel=1, a.wght=5, alpha=1.0)
  lambdaFit<- LKrigFindLambda( x,y,LKinfo=LKinfo)
  lambdaFit$summary


# grid search over parameters 
NG<-3
#NOTE: make this larger ( e.g. 15) for a better (real!) grid search on log lambda
# values
  par.grid<- list( a.wght= rep( 4.05,NG),alpha= rep(1, NG),
                      llambda=  seq(-8,-2,,NG))
  LKinfo<- LKrigSetup(x,NC=12,nlevel=1, a.wght=5, alpha=1.0)
  lambda.search.results<-LKrig.MLE( x,y,LKinfo=LKinfo,
                                    par.grid=par.grid,
                                    lambda.profile=FALSE)
  lambda.search.results$summary
# profile likelihood
  plot( lambda.search.results$summary[,1:2], 
         xlab="effective degrees freedom",
         ylab="ln profile likelihood")
# fit at largest likelihood value:
    lambda.MLE.fit<- LKrig( x,y,
                    LKinfo=lambda.search.results$LKinfo.MLE)
                    
\dontrun{                    
# optimizing  Profile likelihood over lambda using optim
# consider 3 values for a.wght (range parameter)
# in this case the log lambdas passed are the starting values for optim.
  NG<-3
  par.grid<- list( a.wght= c( 4.05,4.1,5),alpha= rep(1, NG),
                      llambda= c(-5,NA,NA))
# NOTE: NAs in llambda mean use the previous MLE for llambda as the
# current starting value. 
  LKinfo<- LKrigSetup(x,NC=12,nlevel=1, a.wght=5, alpha=1.0) 
  lambda.search.results<-LKrig.MLE(
                              x,y,LKinfo=LKinfo, par.grid=par.grid,
                              lambda.profile=TRUE)
  print(lambda.search.results$summary)
# note first result a.wght = 4.05 is the optimized result for the grid
# search given above.
}
########################################################################    
# search over two multi-resolution levels varying the  levels of alpha's
########################################################################
\dontrun{
# NOTE: search ranges found largely by trial and error to make this
# example work also the grid is quite coarse ( and NC is small) to
# be quick as a help file example
  data(NorthAmericanRainfall)
# rename for less typing
  x<- cbind( NorthAmericanRainfall$longitude, NorthAmericanRainfall$latitude)
# total precip in 1/10 mm for JJA 
 y<- log10(NorthAmericanRainfall$precip)
# cut down the size of this data set so examples run quickly
# examples also work with  the full data set. Also try NC= 100 for a
# nontrivial model.
  ind<- x[,1] > -90 & x[,2] < 35 #
  x<- x[ind,]
  y<- y[ind]
  
  Ndes<- 10  
# NOTE: this is set is  very small just to make example run fast
  set.seed(124)
  par.grid<- list()
# create grid of alphas to sum to 1 use a mixture model parametrization
#  alpha1 = (1/(1 + exp(gamma1)) ,
# alpha2 = exp( gamma1) / ( 1 + exp( gamma1))
# 
  par.grid$gamma<- cbind(runif( Ndes, -3,2), runif( Ndes, -3,2))
  par.grid$a.wght<- matrix( 4.5, nrow=Ndes, ncol=3)
# log lambda grid search values
  par.grid$llambda<- runif( Ndes,-5,-3)  
   LKinfo1<- LKrigSetup( x, NC=5, nlevel=3, a.wght=5, alpha=c(1.0,.5,.25))
# NOTE: a.wght in call is not used. Also a better search is to profile over
#  llambda

 alpha.search.results<- LKrig.MLE( x,y,LKinfo=LKinfo1, par.grid=par.grid,
                                    lambda.profile=FALSE)

########################################################################
# Viewing the search results
########################################################################

# this scatterplot is good for a quick look because  effective degrees
# of freedom is a useful summary of fit. 
  plot( alpha.search.results$summary[,1:2], 
         xlab="effective degrees freedom",
         ylab="ln profile likelihood")
#
}

\dontrun{
# a two level model search 
# with profiling over lambda.
#  This takes a few minutes
  Ndes<- 40 
  nlevel<-2 
  par.grid<- list()
## create grid of alphas to sum to 1 use a mixture model parametrization:
#    alpha1 = (1/(1 + exp(gamma1)) ,
#   alpha2 = exp( gamma1) / ( 1 + exp( gamma1))
  set.seed(123)
  par.grid$gamma<- runif( Ndes,-3,4)
## values for range (a.wght)
  a.wght<-  4 + 1/exp(seq( 0,4,,Ndes))
  par.grid$a.wght<- cbind( a.wght, a.wght)
# log lambda grid search values (these are the starting values)
  par.grid$llambda<- rep(-4, Ndes)

  LKinfo1<- LKrigSetup( x, NC=15, nlevel=nlevel, 
                          a.wght=5, alpha=rep( NA,2) ) 
##
## the search over the parameter list in par.grid  maximizing over lambda 
  search.results<- LKrig.MLE( x,y,LKinfo=LKinfo1, par.grid=par.grid,
                                 lambda.profile=TRUE)
# plotting results of likelihood search
set.panel(1,2)
 plot( search.results$summary[,1:2], 
         xlab="effective degrees freedom",
         ylab="ln profile likelihood")
 xtemp<- matrix(NA, ncol=2, nrow=Ndes)
 for( k in 1:Ndes){
   xtemp[k,] <- c( (search.results$par.grid$alpha[[k]])[1],
                  (search.results$par.grid$a.wght[[k]])[1] )
}
 quilt.plot( xtemp,search.results$summary[,2])
#  fit using Tps
 tps.out<- Tps(  xtemp,search.results$summary[,2], lambda=0)
 contour( predictSurface(tps.out), lwd=3,add=TRUE)
 set.panel()
}
\dontrun{
# searching over nu 
data(ozone2)
x<- ozone2$lon.lat
y<- ozone2$y[16,]
good<- !is.na(y)
y<- y[good]
x<- x[good,]
par.grid<- expand.grid( nu=seq(.5,1.5,,4), a.wght=c(4.01,4.1, 4.2,4.5,5))
par.grid$llambda<- rep( NA, length(par.grid$nu))
LKinfo<- LKrigSetup(x,  nlevel=3,NC=5)
out<- LKrig.MLE( x,y, LKinfo=LKinfo, par.grid=par.grid)
}
\dontrun{
# an MLE fit taking advantage of replicated fields

N<- 200
M<-10
sigma<- sqrt(.01)
set.seed(123)
x<- matrix( runif(N*2), N,2)
                
LKinfo<- LKrigSetup( x, NC=16, nlevel=1,
                 a.wght=4.5, lambda=.01, 
                fixed.Function=NULL,
                normalize=TRUE)  
                
# the replicate fields
truef<-  LKrig.sim(x,LKinfo=LKinfo, M=M )
set.seed(222)
error<- sigma*matrix( rnorm(N*M), N,M)
y<- truef + error 
# with correct lambda
obj<- LKrig( x,y, LKinfo=LKinfo, lambda=.01, )
print( obj$sigma.MLE.FULL)
print( obj$rho.MLE.FULL)

fitMLE1<- LKrigFindLambda( x,y, LKinfo=LKinfo)
fitMLE1$summary
aWghtGrid<-  c( 4.01, 4.05, 4.1, 4.2, 4.5, 4.7, 5, 8, 16)
par.grid<- list( a.wght = aWghtGrid)

fitMLE2<- LKrig.MLE( x,y, LKinfo=LKinfo,
                      par.grid= par.grid )
fitMLE2$summary    
plot(  log( aWghtGrid -4)/2, fitMLE2$summary[,2])

LKinfo1<- LKinfoUpdate( LKinfo, lambda=.1, a.wght= 4.2)                   
fitMLE4<- LKrigFindLambdaAwght( x,y, LKinfo=LKinfo1)
fitMLE4$summary[2]

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

