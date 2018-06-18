pkgname <- "LatticeKrig"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('LatticeKrig')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("LKinfo")
### * LKinfo

flush(stderr()); flush(stdout())

### Name: LKinfo
### Title: Specifying the Lattice Krig covariance object('LKinfo') and
###   related utility functions.
### Aliases: LKinfo LKrig.setup LKrig.make.centers LKinfoUpdate
###   LKrig.make.Normalization LKrig.make.a.wght LKrig.make.alpha
###   LKrig.make.grid.info LKrigMRFDecomposition
### Keywords: spatial

### ** Examples

# Load ozone data set
  data(ozone2)  
  x<-ozone2$lon.lat
  y<- ozone2$y[16,]
# Find location that are not 'NA'.
# (LKrig is not set up to handle missing observations.)
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
  LKinfo<- LKrig.setup( x,NC=10,nlevel=2, alpha=c(1,.5),
                            lambda= NA , a.wght=c(5,5))



cleanEx()
nameEx("LKrig.MLE")
### * LKrig.MLE

flush(stderr()); flush(stdout())

### Name: LKrig.MLE
### Title: Simple function to search over covariance parameters for Lattice
###   Krig
### Aliases: LKrig.MLE LKrigFindLambda LKrig.make.par.grid
### Keywords: spatial

### ** Examples

# 

# fitting summer precip for  sub region of North America
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

# This is a single level smoother
  LKinfo<- LKrig.setup(x,NC=12,nlevel=1, a.wght=4.05, alpha=1.0)
  lambdaFit<- LKrigFindLambda( x,y,LKinfo=LKinfo)
  lambdaFit$summary

  NG<-5
 #NOTE: make this larger ( e.g. 15) for a better grid search on log lambda
  par.grid<- list( a.wght= rep( 4.05,NG),alpha= rep(1, NG),
                      llambda=  seq(-8,-2,,NG))
  LKinfo<- LKrig.setup(x,NC=12,nlevel=1, a.wght=5, alpha=1.0)
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
# optimizing  Profile likelihood over lambda using optim
# consider 3 values for a.wght (range parameter)
# in this case the log lambdas passed are the starting values for optim.
  NG<-3
  par.grid<- list( a.wght= c( 4.05,4.1,5),alpha= rep(1, NG),
                      llambda= c(-5,NA,NA))
# NOTE: NAs in llambda mean use the previous MLE for llambda as the
# current starting value. 
  LKinfo<- LKrig.setup(x,NC=12,nlevel=1, a.wght=5, alpha=1.0) 
  lambda.search.results<-LKrig.MLE(
                              x,y,LKinfo=LKinfo, par.grid=par.grid,
                              lambda.profile=TRUE, verbose=TRUE)
  print(lambda.search.results$summary)
# note first result a.wght = 4.05 is the optimized result for the grid
# search given above.

########################################################################    
# search over two multi-resolution levels varying the  levels of alpha's
########################################################################

# NOTE: search ranges found largely by trial and error to make this
# example work also the grid is quite coarse ( and NC is small) to
# be quick as a help file example

  Ndes<- 10  # NOTE: this is set very small just to make example run fast
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
   LKinfo1<- LKrig.setup( x, NC=5, nlevel=3, a.wght=5, alpha=c(1.0,.5,.25))
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
## Not run: 
##D # a more defensible two level model search 
##D # with profiling over lambda.
##D #  This takes a few minutes
##D   Ndes<- 40 
##D   nlevel<-2 
##D   par.grid<- list()
##D ## create grid of alphas to sum to 1 use a mixture model parametrization:
##D #    alpha1 = (1/(1 + exp(gamma1)) ,
##D #   alpha2 = exp( gamma1) / ( 1 + exp( gamma1))
##D   set.seed(123)
##D   par.grid$gamma<- runif( Ndes,-3,4)
##D ## values for range (a.wght)
##D   a.wght<-  4 + 1/exp(seq( 0,4,,Ndes))
##D   par.grid$a.wght<- cbind( a.wght, a.wght)
##D # log lambda grid search values (these are the starting values)
##D   par.grid$llambda<- rep(-4, Ndes)
##D   LKinfo1<- LKrig.setup( x, NC=15, nlevel=nlevel, 
##D                           a.wght=5, alpha=c(1.0,.5,.25) )
##D ##
##D ## the search over the parameter list in par.grid  maximizing over lambda 
##D   search.results<- LKrig.MLE( x,y,LKinfo=LKinfo1, par.grid=par.grid,
##D                                  lambda.profile=TRUE)
##D # plotting results
##D set.panel(1,2)
##D  plot( search.results$summary[,1:2], 
##D          xlab="effective degrees freedom",
##D          ylab="ln profile likelihood")
##D  xtemp<- matrix(NA, ncol=2, nrow=Ndes)
##D  for( k in 1:Ndes){
##D    xtemp[k,] <- c( (search.results$par.grid$alpha[[k]])[1],
##D                   (search.results$par.grid$a.wght[[k]])[1] )
##D }
##D  quilt.plot( xtemp,search.results$summary[,2])
##D #  fit using Tps
##D  tps.out<- Tps(  xtemp,search.results$summary[,2], lambda=0)
##D  contour( predictSurface(tps.out), lwd=3,add=TRUE)
##D  set.panel()
## End(Not run)
## Not run: 
##D # searching over nu 
##D data(ozone2)
##D x<- ozone2$lon.lat
##D y<- ozone2$y[16,]
##D good<- !is.na(y)
##D y<- y[good]
##D x<- x[good,]
##D par.grid<- expand.grid( nu=seq(.5,1.5,,4), a.wght=c(4.01,4.1, 4.2,4.5,5))
##D par.grid$llambda<- rep( NA, length(par.grid$nu))
##D LKinfo<- LKrig.setup(x,  nlevel=3,NC=5)
##D out<- LKrig.MLE( x,y, LKinfo=LKinfo, par.grid=par.grid, verbose=TRUE)
## End(Not run)




cleanEx()
nameEx("LKrig")
### * LKrig

flush(stderr()); flush(stdout())

### Name: LKrig
### Title: Spatial prediction and inference using a compactly supported
###   multi-resolution basis and a lattice model for the basis
###   coefficients.
### Aliases: LKrig predict.LKrig predictSE.LKrig print.LKrig print.LKinfo
###   surface.LKrig predictSurface.LKrig
### Keywords: spatial

### ** Examples

# Load ozone data set
  data(ozone2)  
  x<-ozone2$lon.lat
  y<- ozone2$y[16,]
# Find location that are not 'NA'.
# (LKrig is not set up to handle missing observations.)
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]

# fairly arbitrary choices for covariance parameters and lambda
# just to show a basic level call
  obj1<- LKrig( x,y, a.wght=5, nlevel=3, nu=1.0, NC=4, lambda=.1)

# thin plate spline-like model with the lambda parameter estimated by
# maximum likelihood. Default choices are made for a.wght, nlevel, NC
# and alpha.
## Not run: 
##D   obj<- LatticeKrig( x, y)
##D # summary of fit and a plot of fitted surface
##D   print( obj)
##D   surface( obj )
##D   US(add=TRUE)
##D   points(x)
##D # prediction standard errors
##D   out.se<- predictSE( obj, xnew= x)
## End(Not run)


# breaking down the LatticeKrig function into several steps. 
# also use a different covariance model that has fewer basis fucntions
# (to make the example run more quickly)
 
  LKinfo<- LKrig.setup( x, nlevel=1, alpha=1, NC=15, a.wght=5,
                        lambda=.01)
# maximize likelihood over lambda see help( LKrig.MLE) for details
# this assumes the value of 5 for a.wght.  In many cases the fit is not
# very sensitive to the range parameter such as a.wght in this case --
# but very sensitive to lambda when varied on a log scale.

  MLE.fit<- LKrig.MLE(x,y, LKinfo=LKinfo)
  MLE.fit$summary # summary of optimization over lambda.
# fit using MLE for lambda MLE function has added MLE value of lambda to
# the LKinfo object.
  obj<- LKrig( x,y, LKinfo=MLE.fit$LKinfo.MLE)  
  print( obj)  

# find prediction standard errors at locations based on fixing covariance
# at MLE's
  out.se<- predictSE( obj, xnew= x)
# one could evalute the SE on a grid to get the surface of predicted SE's 
# for large grids it is better to use LKrig.sim.conditional to estimate
#  the variances by Monte Carlo

##########################################################################
# Use multiresolution model that approximates an exponential covariance
# Note that a.wght realted to a range/scale parameter
# is specified at a (seemingly) arbitrary value. 
##########################################################################
  
  LKinfo<- LKrig.setup( x, NC=6, nu=1, nlevel=3, a.wght= 5)
# take a look at the implied covariance function solid=along x
#  and  dashed=along y 
  check.cov<- LKrig.cov.plot( LKinfo)
  matplot( check.cov$d, check.cov$cov, type="l", lty=c(1,2))  

############################################################################
# Search over lambda to find MLE for ozone data with approximate exponential
# covariance
###########################################################################
## Not run: 
##D   
##D   LKinfo.temp<-  LKrig.setup( x, NC=6, nu=1,  nlevel=3, a.wght= 5)
##D # starting value for lambda optimzation
##D   LKinfo.temp$lambda<- 1 
##D   MLE.search<- LKrig.MLE(x,y, LKinfo=LKinfo.temp)
##D # this function returns an LKinfo object with the MLE for lambda included.
##D   MLE.ozone.fit<- LKrig( x,y,  LKinfo= MLE.search$LKinfo.MLE)
## End(Not run) 
###########################################################################
# Including a covariate (linear fixed part in spatial model)
##########################################################################

  data(COmonthlyMet)
  y.CO<- CO.tmin.MAM.climate
  good<- !is.na( y.CO)
  y.CO<-y.CO[good]
  x.CO<- as.matrix(CO.loc[good,])
  Z.CO<- CO.elev[good]
# single level with large range parameter -- similar to a thin plate spline
#  lambda specified 

# fit with elevation
  obj.CO.elev<- LKrig( x.CO,y.CO,Z=Z.CO, nlevel=1, NC=15, alpha=1, lambda=.005,
                          a.wght=4.1)
# BTW the coefficient for the linear term for elevation  is obj.CO$d.coef[4]
# fitted surface without the elevation term
## Not run: 
##D    LKinfo<- LKrig.setup( x.CO, nlevel=1, NC=20,alpha=1, a.wght=4.1, lambda=1.0)
##D   # lambda is the starting vlaue for MLE optimization
##D   CO.MLE<- LKrig.MLE( x.CO,y.CO,Z=Z.CO, LKinfo=LKinfo)
##D   obj.CO.elev<- LKrig( x.CO,y.CO,Z=Z.CO, LKinfo= CO.MLE$LKinfo.MLE)
##D   CO.surface2<- predictSurface( obj.CO.elev, drop.Z=TRUE, nx=50, ny=50)
##D # pull off CO elevations at same locations on grid as the surface
##D   data( RMelevation) # a superset of elevations at 4km resolution
##D   elev.surface<- interp.surface.grid( RMelevation, CO.surface2)
##D    CO.full<- predictSurface( obj.CO.elev, ZGrid= elev.surface, nx=50, ny=50)
##D    
##D # for comparison fit without elevation as a linear covariate:
##D   CO.MLE2<- LKrig.MLE( x.CO,y.CO, LKinfo=LKinfo)
##D   obj.CO<- LKrig( x.CO,y.CO, LKinfo= CO.MLE2$LKinfo.MLE)
##D # surface estimate
##D   CO.surface<- predictSurface( obj.CO, nx=50, ny=50)
##D   set.panel( 2,1)
##D   coltab<- topo.colors(256)
##D   zr<- range( c( CO.full$z), na.rm=TRUE)
##D   image.plot( CO.surface,  col=coltab, zlim =zr)
##D     US( add=TRUE,lwd=2)
##D     title( "MAM min temperatures without elevation")
##D   image.plot( CO.full, col=coltab, zlim=zr)
##D     title( "Including elevation")
##D     US( add=TRUE,lwd=2)
## End(Not run)

########################################################################
# for a more elaborate search over  a.wght, alpha and lambda to find
# joint MLEs see help(LKrig.MLE)
########################################################################

########################################################################
# A bigger problem: 26K observations and 4.6K basis functions
# fitting takes about 15 seconds on a laptop for a fixed covariance
#  LKrig.MLE to find the MLE (not included) for lambda takes abou
#  8 minutes
#######################################################################
## Not run: 
##D   data(CO2)
##D   obj1<- LKrig( CO2$lon.lat,CO2$y,NC=100,nlevel=1, lambda=.088,
##D                        a.wght=5, alpha=1)
##D # 4600 basis functions 100X46 lattice  (number of basis functions
##D # reduced in y direction because of a rectangular domain
##D   obj1$trA.est # about 1040 effective degrees of freedom 
##D #
##D   glist<- list( x= seq( -180,180,,200),y=seq( -80,80,,100) )
##D   xg<-  make.surface.grid(glist)
##D   fhat<- predict( obj1, xg)
##D   fhat <- matrix( fhat,200,100) # convert to image
##D #Plot data and gap-filled estimate
##D   set.panel(2,1)
##D   quilt.plot(CO2$lon.lat,CO2$y,zlim=c(373,381))
##D   title("Simulated CO2 satellite observations")
##D   world(add=TRUE,col="magenta")
##D   image.plot( glist$x, glist$y, fhat,zlim=c(373,381))
##D   world( add=TRUE, col="magenta")
##D   title("Gap-filled global predictions")
## End(Not run)  
## same example but use a periodic lattice in longitude and 
## periodic basis functions. This is the "cylinder" model.
## Not run: 
##D   data(CO2)
##D   LKinfo.cyl<- LKrig.setup( cbind( c( -180,180), range( CO2$lon.lat[,2])),
##D                             NC=25,nlevel=1, lambda=.1,
##D                              a.wght=5, alpha=1,  distance.type="cylinder")
##D # 1127 basis functions from a 49X23 lattice
##D   search.CO2<- LKrig.MLE( CO2$lon.lat,CO2$y, LKinfo=LKinfo.cyl)
##D # MLE search over lambda
##D   obj2<- LKrig( CO2$lon.lat,CO2$y, LKinfo=search.CO2$LKinfo.MLE)
##D   surface( obj2)
##D   world( add=TRUE)
## End(Not run)  

set.panel()
#########################################################################
#  Comparing LKrig to ordinary Kriging
########################################################################

# Here is an illustration of using the fields function mKrig with the
# LKrig covariance to reproduce the computations of LKrig. The
# difference is that mKrig can not take advantage of any sparsity in
# the precision matrix because its inverse, the covariance matrix, is
# not sparse.  This example reinforces the concept that LKrig finds the
# the standard geostatistical estimate but just uses a particular
# covariance function defined via basis functions and the precision
# matrix.
# Load ozone data set (AGAIN)
## Not run: 
##D   data(ozone2)  
##D   x<-ozone2$lon.lat
##D   y<- ozone2$y[16,]
##D # Find location that are not 'NA'.
##D # (LKrig is not set up to handle missing observations.)
##D   good <-  !is.na( y)
##D   x<- x[good,]
##D   y<- y[good]
##D   a.wght<- 5
##D   lambda <-  1.5
##D   obj1<- LKrig( x,y,NC=16,nlevel=1, alpha=1,  lambda=lambda, a.wght=5,
##D                 NtrA=20,iseed=122)
##D  
##D # in both calls iseed is set the same so we can compare 
##D # Monte Carlo estimates of effective degrees of freedom
##D   obj1$trA.est
##D # The covariance "parameters" are all in the list LKinfo
##D # to create this special list outside of a call to LKrig use
##D   LKinfo.test <- LKrig.setup( x, NC=16, nlevel=1, alpha=1.0,  a.wght=5)
##D 
##D # this call to mKrig should be identical to the LKrig results
##D # because it uses the LKrig.cov covariance with all the right parameters.
##D   obj2<- mKrig( x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
##D                       cov.args=list( LKinfo=LKinfo.test), NtrA=20, iseed=122)
##D # compare the two results this is also an
##D # example of how tests are automated in fields
##D # set flag to have tests print results
##D   test.for.zero.flag<- TRUE
##D   test.for.zero( obj1$fitted.values, obj2$fitted.values,
##D                   tag="comparing predicted values LKrig and mKrig")
##D # compare standard errors. 
##D   se1<- predictSE.LKrig( obj1)
##D   se2<- predictSE.mKrig(obj2)
##D   test.for.zero( se1, se2,
##D                   tag="comparing standard errors for LKrig and mKrig")
## End(Not run)
########################################################################
#  Unnormalized marginal variance for a 6X6 basis on [-1,1]X[-1,1]
#  This is an example of why normalization seems important. 
########################################################################

## Not run: 
##D # specify covariance without normalization note all we need is the
##D #corners of domains to setup the info list. 
##D   LKinfo<- LKrig.setup(cbind( c( -1,1), c(-1,1)), NC=6, nlevel=1,
##D                         a.wght=4.5,alpha=1, normalize= FALSE)  
##D # 80X80 grid of points 
##D   xg<- make.surface.grid( list(x=seq(-1,1,,80), y= seq(-1,1,,80))) 
##D   look<- LKrig.cov( xg, LKinfo=LKinfo,marginal =TRUE) 
##D # surface plot of the marginal variances of the process. 
##D   image.plot( as.surface(xg, look)) 
##D # basis function centers from the first (and only) level
##D   xp<- make.surface.grid( LKinfo$grid[[1]])
##D   points(xp)
## End(Not run)




cleanEx()
nameEx("LKrig.basis")
### * LKrig.basis

flush(stderr()); flush(stdout())

### Name: LKrig.basis
### Title: Functions for generating a multiresolution, compactly supported
###   basis, multiresolution covariance functions and simulating from these
###   processes.
### Aliases: LKrig.basis LKrig.precision LKrig.MRF.precision LKrig.cov
###   LKrig.normalize.basis LKrig.normalize.basis.fast LKrig.cov.plot
###   LKrig.quadraticform LKrig.spind2spam
### Keywords: spatial

### ** Examples

# Load ozone data set
  data(ozone2)  
  x<-ozone2$lon.lat
  y<- ozone2$y[16,]
# Find location that are not 'NA'.
# (LKrig is not set up to handle missing observations.)
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
  LKinfo<- LKrig.setup( x,NC=20,nlevel=1, alpha=1, lambda= .3 , a.wght=5)
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
    LKinfo2<- LKrig.setup( cbind( c( -1,1), c(-1,1)), NC=NC,
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
## Not run: 
##D # full 2-d view of the covariance (this example follows the code
##D # in LKrig.cov.plot)
##D  x2<- cbind( 0,0)
##D  x1<- make.surface.grid( list(x=seq( -1,1,,40),  y=seq( -1,1,,40)))
##D  look<- LKrig.cov( x1,x2, LKinfo2)
##D  contour( as.surface( x1, look))
##D # Note nearly circular contours.
##D # of  course  plot(look[,80/2]) should look like plot above.
##D #
## End(Not run)

## Not run: 
##D #Some correlation functions from different models
##D set.panel(2,1)
##D # a selection of ranges:
##D   hold<- matrix( NA, nrow=150, ncol=4)
##D   kappa<- seq( .25,1,,4)
##D   x2<- cbind( 0,0)
##D   x1<-  cbind( seq(-1,1,,150), rep( 0,150))
##D   for( k in 1:4){
##D     LKtemp<-  LKrig.setup( cbind( c( -1,1), c(-1,1)), NC=NC,
##D                    nlevel=nlevel,
##D                    a.wght= 4  + 1/(kappa[k]^2),
##D                    alpha=alpha)
##D     hold[,k]<-  LKrig.cov( x1,x2, LKinfo=LKtemp)
##D   }
##D   matplot( x1[,1], hold, type="l", lty=1, col=rainbow(5), pch=16 )
##D # a selection of smoothness parameters
##D   ktemp<- .5 # fix range
##D   alpha.power<- seq( 1,4,,4)
##D   LKtemp<- LKinfo2
##D   for( k in 1:4){
##D    LKtemp<-  LKrig.setup( cbind( c( -1,1), c(-1,1)), NC=NC,
##D                    nlevel=nlevel,
##D                    a.wght= 4  + 1/(ktemp^2),
##D                    alpha=alpha^alpha.power[k])
##D     hold[,k]<-  LKrig.cov( x1,x2, LKinfo=LKtemp)
##D   }
##D   matplot( x1[,1], hold, type="l", lty=1, col=rainbow(5) )
##D  set.panel()
## End(Not run)
 
## Not run: 
##D # generating a basis on the domain [-1,1] by [-1,1] with 1 level
##D # Default number of buffer points are added to each side. 
##D   LKinfo<- LKrig.setup(cbind( c(-1,1), c(-1,1)), NC=6,
##D                                  nlevel=1, a.wght=4.5,alpha=1, NC.buffer=0 )
##D # evaluate the basis functions on a grid to look at them
##D   xg<- make.surface.grid( list(x=seq(-1,1,,50), y= seq(-1,1,,50)))
##D   PHI<- LKrig.basis( xg,LKinfo)
##D   dim(PHI) # should be  2500=50^2  by  36=6^2
##D # plot the 9th basis function  as.surface is a handy function to
##D # reformat the vector as an image object
##D # using the grid information in an attribute of the grid points 
##D   image.plot(as.surface(xg, PHI[,9]))
##D   points(  make.surface.grid( LKinfo$grid[[1]]), col="grey", cex=.5)
##D 
##D set.panel()
## End(Not run)
#
# example of basis function indexing
#
## Not run: 
##D # generating a basis on the domain [-1,1]X[-1,1] with 3 levels
##D # note that there are no buffering grid points.
##D   set.panel(3,2)
##D   LKinfo<-LKrig.setup(cbind( c(-1,1), c(-1,1)), NC=6,
##D                     a.wght=rep(5,3), alpha=c(1,.5,.25), nlevel=3,
##D                     NC.buffer=0)
##D # evaluate the basis functions on a grid to look at them
##D   xtemp<- seq(-1,1,,40)
##D   xg<- make.surface.grid( list(x=xtemp, y= xtemp) )
##D   PHI<- LKrig.basis( xg,LKinfo)
##D # coerce to dense matrix format to make plotting easier.
##D   PHI<- spam2full(PHI)
##D # first tenth, and last basis function in each resolution level
##D # basis functions centers are added
##D   set.panel(3,3)
##D   grid.info<- LKinfo$grid.info
##D   for(  j in 1:3){
##D     id1<- LKinfo$offset[j]+ 1
##D     id2<-  LKinfo$offset[j]+ 10
##D     idlast<- LKinfo$offset[j]+ LKinfo$mx[j]*LKinfo$my[j]
##D  
##D     centers<-  make.surface.grid( LKinfo$grid[[j]] )
##D     image.plot( as.surface(xg, PHI[,id1]))
##D     points( centers, cex=.2, col="grey")
##D     image.plot(as.surface(xg, PHI[,id2]))
##D     points( centers, cex=.2, col="grey")
##D     image.plot( as.surface(xg, PHI[,idlast]))
##D     points( centers, cex=.2, col="grey")}
##D 
##D   set.panel()
## End(Not run)
## Not run: 
##D # examining the stationarity of covariance model
##D   temp.fun<- 
##D      function( NC.buffer=0, NC=4,  a.wght=4.01){
##D         LKinfo<- LKrig.setup(cbind( c(-1,1), c(-1,1)),nlevel=1, alpha=1,
##D                                  a.wght=a.wght, NC=NC,  NC.buffer=NC.buffer)
##D         cov1y<- cov1x<- cov0x<- cov0y<-  matrix( NA, nrow=200, ncol=20)
##D         cov1dx<- cov1dy<- cov0dx<- cov0dy<- matrix( NA, nrow=200, ncol=20)
##D         cgrid<- seq( 0,1,,20)
##D         for( k in 1:20){
##D             hold<- LKrig.cov.plot( LKinfo,
##D                             center=rbind( c(cgrid[k], cgrid[k])), NP=200)
##D             cov1x[,k] <- hold$cov[,1]
##D             cov1y[,k] <- hold$cov[,2]
##D             cov1dx[,k] <- hold$d[,1]
##D             cov1dy[,k] <- hold$d[,2]
##D             hold<- LKrig.cov.plot( LKinfo,
##D                              center=rbind( c(cgrid[k],0) ), NP=200)
##D             cov0x[,k] <- hold$cov[,1]
##D             cov0y[,k] <- hold$cov[,2]
##D             cov0dx[,k] <- hold$d[,1]
##D             cov0dy[,k] <- hold$d[,2]
##D                 }
##D          matplot( cov1dx, cov1x, type="l", col= rainbow(20),
##D                          xlab="", ylab="correlation")
##D          mtext( side=1, line=-1, text="diagonal X")
##D          title( paste(  " buffer=",NC.buffer), cex=.5)
##D          matplot( cov1dy, cov1y, type="l", col= rainbow(20),
##D                         xlab="", ylab="")
##D          mtext( side=1, line=-1, text="diagonal Y")
##D          matplot(cov0dx, cov0x, type="l", col= rainbow(20),
##D                         xlab="",       ylab="")
##D          mtext( side=1, line=-1, text="middle X")
##D          matplot( cov0dy, cov0y, type="l", col= rainbow(20),
##D                          xlab="",   ylab="")
##D          mtext( side=1, line=-1, text="middle Y")
##D          title( paste( NC, a.wght), cex=.5)
##D }
##D 
##D 
##D  set.panel(3,4)
##D par(mar=c(3,4,1,0), oma=c(1,1,1,1))
##D temp.fun(  NC.buffer=5, NC=4, a.wght=4.05)
##D temp.fun(  NC.buffer=5, NC=16, a.wght=4.05)
##D temp.fun(  NC.buffer=5, NC=64, a.wght=4.05)
##D 
##D set.panel(4,4)
##D par(mar=c(3,4,1,0), oma=c(1,1,1,1))
##D temp.fun( NC.buffer=0, NC=8)
##D temp.fun( NC.buffer=2, NC=8)
##D temp.fun( NC.buffer=4, NC=8)
##D # this next one takes a while
##D temp.fun( NC.buffer=8,  NC=8)
##D # stationary == curves should all line up!
##D 
## End(Not run)




cleanEx()
nameEx("LKrig.fixed")
### * LKrig.fixed

flush(stderr()); flush(stdout())

### Name: LKrigDefaultFixedFunction
### Title: Creates fixed part of spatial model.
### Aliases: LKrigDefaultFixedFunction predictLKrigFixedFunction
### Keywords: spatial

### ** Examples

x<- matrix( runif(100), nrow=50)
# linear polynomial 
T.matrix<- LKrigDefaultFixedFunction(x, m=2)
# quadratic polynomial 
T.matrix<- LKrigDefaultFixedFunction(x, m=3)



cleanEx()
nameEx("LKrig.sim")
### * LKrig.sim

flush(stderr()); flush(stdout())

### Name: LKrig.sim
### Title: Functions for simulating a multiresolution process following the
###   Lattice Krig covariance model.
### Aliases: LKrig.sim LKrig.sim.conditional
### Keywords: spatial

### ** Examples

# Load ozone data set
  data(ozone2)  
  x<-ozone2$lon.lat
  y<- ozone2$y[16,]
# Find location that are not 'NA'.
# (LKrig is not set up to handle missing observations.)
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
  LKinfo<- LKrig.setup( x,NC=20,nlevel=1, alpha=1, lambda= .3 , a.wght=5)
# BTW lambda is close to MLE 
# Simulating this  LKrig process
# simulate 4 realizations of process and plot them
# (these have unit marginal variance)
  xg<- make.surface.grid(list( x=seq( -87,-83,,40), y=seq(36.5, 44.5,,40)))
  out<- LKrig.sim(xg, LKinfo,M=4)
## Not run: 
##D   set.panel(2,2)
##D   for( k in 1:4){
##D     image.plot( as.surface( xg, out[,k]), axes=FALSE) }
## End(Not run)
  obj<- LKrig(x,y,LKinfo=LKinfo)
  O3.cond.sim<- LKrig.sim.conditional( obj, M=3,nx=40,ny=40) 
## Not run: 
##D   set.panel( 2,2)
##D   zr<- range( c(  O3.cond.sim$draw,  O3.cond.sim$ghat), na.rm=TRUE)
##D   coltab<- tim.colors()
##D   image.plot( as.surface( O3.cond.sim$x.grid, O3.cond.sim$ghat), zlim=zr)
##D   title("Conditional mean")
##D   US( add=TRUE)
##D   for( k in 1:3){
##D     image( as.surface( O3.cond.sim$x.grid, O3.cond.sim$g.draw[,k]),
##D               zlim=zr, col=coltab)
##D     points( obj$x, cex=.5)
##D     US( add=TRUE)
##D   }
##D   set.panel()
## End(Not run)




cleanEx()
nameEx("LatticeKrig")
### * LatticeKrig

flush(stderr()); flush(stdout())

### Name: LatticeKrig
### Title: User-friendly spatial prediction and inference using a compactly
###   supported multi-resolution basis and a lattice model for the basis
###   coefficients.
### Aliases: LatticeKrig
### Keywords: spatial

### ** Examples

# Load ozone data set
  data(ozone2)  
  x<-ozone2$lon.lat
  y<- ozone2$y[16,]

# thin plate spline-like model with the lambda parameter estimated by
# maximum likelihood. Default choices are made for a.wght, nlevel, NC
# and alpha.
  obj<- LatticeKrig( x, y)
## Not run: 
##D # summary of fit and a plot of fitted surface
##D   print( obj)
##D   surface( obj )
##D   US(add=TRUE)
##D   points(x)
##D # prediction standard errors
##D   out.se<- predictSE( obj, xnew= x)
## End(Not run)

###########################################################################
# Including a covariate (linear fixed part in spatial model)
########################################################################## 
## Not run: 
##D   data(COmonthlyMet)
##D 
##D   obj  <- LatticeKrig(CO.loc,  CO.tmin.MAM.climate, Z=CO.elev)
##D   obj2 <- LatticeKrig(CO.loc, CO.tmin.MAM.climate)
##D 
##D # compare with and without linear covariates
##D   set.panel(1,2)
##D   surface(obj)
##D   US(add=TRUE)
##D   title("With Elevation Covariate")
##D 
##D   surface(obj2)
##D   US(add=TRUE)
##D   title("Without Elevation Covariate")
##D 
## End(Not run)
## Not run: 
##D  data(COmonthlyMet)
##D # Examining a few different "range" parameters
##D a.wghtGrid<-  4  +  c( .1,.5, .8, 1, 2)^2
##D 
##D #NOTE smallest is "spline like" the largest is essentially independent
##D # coefficients at each level.  In this case the "independent" end is
##D # favored but the eff df. of the surface is very similar across models
##D # indicating about the same separate of the estimates into spatial
##D # signal and noise
##D #
##D for( k in 1:5 ){
##D obj  <- LatticeKrig(CO.loc,  CO.tmin.MAM.climate, Z=CO.elev, 
##D                       a.wght=a.wghtGrid[k])
##D cat( "a.wght:", a.wghtGrid[k], "ln Profile Like:",
##D             obj$lnProfileLike, "Eff df:", obj$trA.est, fill=TRUE)
##D }
## End(Not run)




cleanEx()
nameEx("Radial.Basis")
### * Radial.Basis

flush(stderr()); flush(stdout())

### Name: Radial.basis
### Title: Two dimensional radial basis functions based on a Wendland
###   function.
### Aliases: Radial.basis LKrig.cyl WendlandFunction
### Keywords: spatial

### ** Examples

x<- cbind( runif(100), runif(100))
center<- expand.grid( seq( 0,1,,5), seq(0,1,,5))
# coerce to matrix
center<- as.matrix(center)

PHI<- Radial.basis(x, center, delta=.5)

# LKrig with a different radial basis function. 
# 
  data(ozone2)  
  x<-ozone2$lon.lat
  y<- ozone2$y[16,]
# Find location that are not 'NA'.
# (LKrig is not set up to handle missing observations.)
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
  obj<- LKrig(x,y,NC=30,nlevel=1, alpha=1, lambda=.01, a.wght=5)
  triweight<- function( d){
       ifelse( abs(d)<=1, (1-d^2)^3, 0)}   
  obj1<- LKrig(x,y,NC=30,nlevel=1, alpha=1, 
    lambda=.01, a.wght=5, RadialBasisFunction="triweight", overlap=1.8)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
