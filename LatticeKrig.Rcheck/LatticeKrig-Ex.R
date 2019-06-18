pkgname <- "LatticeKrig"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "LatticeKrig-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('LatticeKrig')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("IcosohedronGrid")
### * IcosohedronGrid

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: IcosohedronGrid
### Title: Icosoherdal multiresolution grids
### Aliases: IcosohedronGrid IcosohedronFaces
### Keywords: spatial

### ** Examples

# second level in lon lat coordinates 
look<- IcosohedronGrid(3)
lonlat<- toSphere( look[[3]])
plot( lonlat, xlab="lon", ylab="lat")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("IcosohedronGrid", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LKDist")
### * LKDist

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LKDist
### Title: Find all pairwise distances within a maximum distance.
### Aliases: LKDist LKDistComponents LKDistGrid LKDistGridComponents
###   LKGridCheck LKGridFindNmax
### Keywords: spatial

### ** Examples

set.seed( 123)	
x<- matrix( runif(100*2), 100,2)

DMatrix<- LKDist( x,x, delta=.1)
# coerce to spam matrix format
DMatrix2<- spind2spam( DMatrix)

# a grid
gridL<- list( x1= seq(0,1,.2), x2= seq( 0,2,.2) , x3= seq( -1,1,.2))
class(gridL)<- "gridList"	
x1<- cbind( runif( 100), runif(100)*2, 2*(runif( 100) -.5) )
look<- LKDistGrid( x1, gridL, delta=.45)
# check against rdist.
# look2<- rdist( x1, make.surface.grid(gridL))
# look2[ look2 >= .45] <- 0
# max( abs(look- look2)[look>0] )

# test of periodic option
 gridL<- structure(
            list( x1= seq(0,1,.02),
                  x2= seq( 0,1,.02)),
            class="gridList")
 look1<- LKDistGrid( rbind(c(0,0)), gridL, delta=.35,
                     periodic=c(TRUE,FALSE))
 look2<- spind2full(look1)
 image.plot( as.surface( gridL, look2) )
 
 look1<- LKDistGrid( rbind(c(0,0)), gridL, delta=.35,
                      periodic=c(TRUE,TRUE))
 look2<- spind2full(look1)
 image.plot( as.surface( gridL, look2) )  



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LKDist", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LKRectangle")
### * LKRectangle

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LKRectangle
### Title: Summary of the LKRectangle geometry for a standard two
###   dimensional spatial domain.
### Aliases: LKRectangle
### Keywords: spatial

### ** Examples

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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LKRectangle", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LKrig.MLE")
### * LKrig.MLE

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LKrig.MLE
### Title: Simple function to search over covariance parameters for Lattice
###   Krig.
### Aliases: LKrig.MLE LKrigFindLambda LKrigFindLambdaAwght
###   LambdaAwghtObjectiveFunction LKrig.make.par.grid omega2Awght
###   Awght2Omega
### Keywords: spatial

### ** Examples

# 
# fitting summer precip for  sub region of North America (Florida)
# (tiny subregion is just to make this run under 5 seconds). 
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
 
  LKinfo<- LKrigSetup(x,NC=4, nlevel=1, a.wght=5, alpha=1.0)
  lambdaFit<- LKrigFindLambda( x,y,LKinfo=LKinfo)
  lambdaFit$summary

## Not run: 
##D # grid search over parameters 
##D   NG<-15
##D   par.grid<- list( a.wght= rep( 4.05,NG),alpha= rep(1, NG),
##D                       llambda=  seq(-8,-2,,NG))
##D   lambda.search.results<-LKrig.MLE( x,y,LKinfo=LKinfo,
##D                                     par.grid=par.grid,
##D                                     lambda.profile=FALSE)
##D   lambda.search.results$summary
##D # profile likelihood
##D   plot( lambda.search.results$summary[,1:2], 
##D          xlab="effective degrees freedom",
##D          ylab="ln profile likelihood")
##D # fit at largest likelihood value:
##D   lambda.MLE.fit<- LKrig( x,y,
##D                     LKinfo=lambda.search.results$LKinfo.MLE)
## End(Not run)                    
                    
## Not run: 
##D                     
##D # optimizing  Profile likelihood over lambda using optim
##D # consider 3 values for a.wght (range parameter)
##D # in this case the log lambdas passed are the starting values for optim.
##D   NG<-3
##D   par.grid<- list( a.wght= c( 4.05,4.1,5) ,alpha= rep(1, NG),
##D                       llambda= c(-5,NA,NA))
##D # NOTE: NAs in llambda mean use the previous MLE for llambda as the
##D # current starting value. 
##D   LKinfo<- LKrigSetup(x,NC=12,nlevel=1, a.wght=5, alpha=1.0) 
##D   lambda.search.results<-LKrig.MLE(
##D                               x,y,LKinfo=LKinfo, par.grid=par.grid,
##D                               lambda.profile=TRUE)
##D   print(lambda.search.results$summary)
##D # note first result a.wght = 4.05 is the optimized result for the grid
##D # search given above.
## End(Not run)
########################################################################    
# search over two multi-resolution levels varying the  levels of alpha's
########################################################################
## Not run: 
##D # NOTE: search ranges found largely by trial and error to make this
##D # example work also the grid is quite coarse ( and NC is small) to
##D # be quick as a help file example
##D   data(NorthAmericanRainfall)
##D # rename for less typing
##D   x<- cbind( NorthAmericanRainfall$longitude, NorthAmericanRainfall$latitude)
##D # total precip in 1/10 mm for JJA 
##D  y<- log10(NorthAmericanRainfall$precip)
##D # cut down the size of this data set so examples run quickly
##D # examples also work with  the full data set. Also try NC= 100 for a
##D # nontrivial model.
##D   ind<- x[,1] > -90 & x[,2] < 35 #
##D   x<- x[ind,]
##D   y<- y[ind]
##D   
##D   Ndes<- 10  
##D # NOTE: this is set to be very small just to make this
##D #       example run fast
##D   set.seed(124)
##D   par.grid<- list()
##D # create grid of alphas to sum to 1 use a mixture model parametrization
##D #  alpha1 = (1/(1 + exp(gamma1)) ,
##D # alpha2 = exp( gamma1) / ( 1 + exp( gamma1))
##D # 
##D   par.grid$gamma<- cbind(runif( Ndes, -3,2), runif( Ndes, -3,2))
##D   par.grid$a.wght<- rep( 4.5, Ndes)
##D # log lambda grid search values
##D   par.grid$llambda<- runif( Ndes,-5,-3)  
##D   LKinfo1<- LKrigSetup( x, NC=5, nlevel=3, a.wght=5, alpha=c(1.0,.5,.25))
##D # NOTE: a.wght in call is not used. Also a better search is to profile over
##D #  llambda
##D 
##D  alpha.search.results<- LKrig.MLE( x,y,LKinfo=LKinfo1, par.grid=par.grid,
##D                                     lambda.profile=FALSE)
##D 
##D ########################################################################
##D # Viewing the search results
##D ########################################################################
##D 
##D # this scatterplot is good for a quick look because  effective degrees
##D # of freedom is a useful summary of fit. 
##D   plot( alpha.search.results$summary[,1:2], 
##D          xlab="effective degrees freedom",
##D          ylab="ln profile likelihood")
##D #
## End(Not run)

## Not run: 
##D # a two level model search 
##D # with profiling over lambda.
##D data(NorthAmericanRainfall)
##D # rename for less typing
##D   x<- cbind( NorthAmericanRainfall$longitude,
##D              NorthAmericanRainfall$latitude)
##D # mean total precip in 1/10 mm for JJA 
##D   y<- log10(NorthAmericanRainfall$precip)
##D 
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
##D   par.grid$a.wght<- 4 + 1/exp(seq( 0,4,,Ndes))
##D # log lambda grid search values (these are the starting values)
##D   par.grid$llambda<- rep(-4, Ndes)
##D 
##D   LKinfo1<- LKrigSetup( x, NC=15, nlevel=nlevel, 
##D                           a.wght=5, alpha=rep( NA,2) ) 
##D ##
##D ## the search over the parameter list in par.grid  maximizing over lambda 
##D   search.results<- LKrig.MLE( x,y,LKinfo=LKinfo1, par.grid=par.grid,
##D                                  lambda.profile=TRUE)
##D # plotting results of likelihood search
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
##D par.grid<- expand.grid( nu = c(.5,1.0, 1.5), a.wght= list(4.1,4.5,5) )
##D par.grid$llambda<- rep( NA, length(par.grid$nu))
##D LKinfo<- LKrigSetup(x,  nlevel=3, nu=.5, NC=5, a.wght=4.5)
##D out<- LKrig.MLE( x,y, LKinfo=LKinfo, par.grid=par.grid)
##D # take a look
##D cbind( par.grid, out$summary[,1:5])
## End(Not run)
## Not run: 
##D # an MLE fit taking advantage of replicated fields
##D # check based on simulated data
##D 
##D N<- 200
##D M<-50 # number of independent replicated fields
##D sigma<- sqrt(.01)
##D set.seed(123)
##D x<- matrix( runif(N*2), N,2)
##D                 
##D LKinfo<- LKrigSetup( x, NC=16, nlevel=1,
##D                  a.wght=4.5, lambda=.01, 
##D                 fixed.Function=NULL,
##D                 normalize=TRUE)  
##D                 
##D # the replicate fields
##D truef<-  LKrig.sim(x,LKinfo=LKinfo, M=M )
##D set.seed(222)
##D error<- sigma*matrix( rnorm(N*M), N,M)
##D y<- truef + error 
##D # with correct lambda
##D obj<- LKrig( x,y, LKinfo=LKinfo, lambda=.01, )
##D print( obj$sigma.MLE.FULL)
##D print( obj$rho.MLE.FULL)
##D 
##D fitMLE1<- LKrigFindLambda( x,y, LKinfo=LKinfo)
##D fitMLE1$summary
##D aWghtGrid<-  c( 4.01, 4.05, 4.1, 4.2, 4.5, 4.6, 4.7, 5, 8, 16)
##D par.grid<- list( a.wght = aWghtGrid)
##D 
##D fitMLE2<- LKrig.MLE( x,y, LKinfo=LKinfo,
##D                       par.grid= par.grid )
##D fitMLE2$summary   
##D 
##D LKinfo1<- LKinfoUpdate( LKinfo, lambda=.1, a.wght= 4.2)                   
##D fitMLE4<- LKrigFindLambdaAwght( x,y, LKinfo=LKinfo1)
##D fitMLE4$summary
##D 
##D plot(  log( aWghtGrid -4)/2, fitMLE2$summary[,2], type="b",
##D   xlab="log( a.wght - 4)/2",
##D   ylab= "log Profile likelihood" )
##D 
##D 
##D points( log(fitMLE4$summary["a.wght.MLE"] -4)/2,
##D      fitMLE4$summary["lnProfLike"], pch="+", col="red"  )
##D xline( log(fitMLE4$summary["a.wght.MLE"] -4)/2, col="red", lty=2)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LKrig.MLE", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LKrig")
### * LKrig

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LKrig
### Title: Spatial prediction and inference using a compactly supported
###   multi-resolution basis and a lattice model for the basis
###   coefficients.
### Aliases: LKrig choleskyMemory predict.LKrig predictSE.LKrig print.LKrig
###   summary.LKrig print.LKinfo surface.LKrig predictSurface.LKrig
### Keywords: spatial

### ** Examples

# NOTE: CRAN requires that the "run" examples  execute within  5 seconds. 
# so to meet this constraint many of these have been 
# severely limited to be quick, typically making NC and nlevel
# very small.
	
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
  obj1<- LKrig( x,y, a.wght=5, nlevel=3, nu=1.0, NC=10, lambda=.1)
  print(obj1)
#  
# this is the same as:
#  LKinfoEX<- LKrigSetup( x, a.wght=5, nlevel=3, nu=1.0, NC=4)
#  obj1<- LKrig( x,y, LKinfo= LKinfoEX, lambda=.1)  

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
# also use a different covariance model that has fewer basis functions
# (to make the example run more quickly)

## Not run: 
##D  
##D   LKinfo<- LKrigSetup( x, nlevel=3, nu=1, NC=5, a.wght=5,
##D                         lambda=.01)
##D # maximize likelihood over lambda see help( LKrig.MLE) for details
##D # this assumes the value of 5 for a.wght.  In many cases the fit is not
##D # very sensitive to the range parameter such as a.wght in this case --
##D # but very sensitive to lambda when varied on a log scale.
##D 
##D   MLE.fit<- LKrig.MLE(x,y, LKinfo=LKinfo)
##D   MLE.fit$summary # summary of optimization over lambda.
##D # fit using MLE for lambda MLE function has added MLE value of lambda to
##D # the LKinfo object.
##D   obj<- LKrig( x,y, LKinfo=MLE.fit$LKinfo.MLE)  
##D   print( obj)  
##D # find prediction standard errors at locations based on fixing covariance
##D # at MLE's
##D   out.se<- predictSE( obj, xnew= x)
##D # one could evaluate the SE on a grid to get the surface of predicted SE's 
##D # for large grids it is better to use LKrig.sim.conditional to estimate
##D #  the variances by Monte Carlo
## End(Not run)

##########################################################################
# Use multiresolution model that approximates an exponential covariance
# Note that a.wght, related to a range/scale parameter
# is specified at a (seemingly) arbitrary value. 
##########################################################################
  
  LKinfo<- LKrigSetup( x, NC=4, nu=1, nlevel=2, a.wght= 5)
# take a look at the implied covariance function solid=along x
#  and  dashed=along y 
  check.cov<- LKrig.cov.plot( LKinfo)
  matplot( check.cov$d, check.cov$cov, type="l", lty=c(1,2))  

############################################################################
# Search over lambda to find MLE for ozone data with approximate exponential
# covariance
###########################################################################
## Not run: 
##D   LKinfo.temp<-  LKrigSetup( x, NC=6, nu=1,  nlevel=3, a.wght= 5, lambda=.5)
##D # starting value for lambda optimization 
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
# note how for convenience the LKrig parameters are
# just included here and not passed as a separate LKinfo object. 
  obj.CO.elev<- LKrig( x.CO,y.CO,Z=Z.CO, nlevel=1, NC=50, alpha=1, lambda=.005,
                          a.wght=4.1)
# BTW the coefficient for the linear term for elevation  is obj.CO$d.coef[4]
# fitted surface without the elevation term
## Not run: 
##D    LKinfo<- LKrigSetup( x.CO, nlevel=1, NC=20,alpha=1, a.wght=4.1, lambda=1.0)
##D # lambda given here is just the starting value for MLE optimization
##D   CO.MLE<- LKrig.MLE( x.CO,y.CO,Z=Z.CO, LKinfo=LKinfo)
##D   obj.CO.elev<- LKrig( x.CO,y.CO,Z=Z.CO, LKinfo= CO.MLE$LKinfo.MLE)
##D   CO.surface2<- predictSurface( obj.CO.elev, drop.Z=TRUE, nx=50, ny=50)
##D # pull off CO elevations at same locations on grid as the surface
##D   data( RMelevation) 
##D # a superset of elevations at 4km resolution
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
## Not run: 
##D #### a 3-d example using alternative geometry
##D set.seed( 123)
##D N<- 500
##D x<-  matrix( runif(3* N,-1,1), ncol=3, nrow=N)
##D y<-   10*exp( -rdist( x, rbind( c(.5,.5,.6) ) )/.5)
##D LKinfo<- LKrigSetup( x,  
##D                  LKGeometry = "LKBox",
##D                      nlevel = 1,
##D                      a.wght = 6.01,
##D                          NC = 5,
##D                   NC.buffer = 2,
##D                   normalize = TRUE,
##D              choleskyMemory = list(nnzR= 2e6),
##D                     lambda = .1,
##D                     mean.neighbor=200
##D                    )  
##D # NOTE: memory for the spam routines has been increased to 
##D # avoid warnings  
##D   out1<- LKrig( x,y, LKinfo=LKinfo)
## End(Not run)  
  
## Not run: 
##D # MLE search over lambda and  a bigger problem
##D # but no normalization
##D N<- 5e4
##D x<-  matrix( runif(3* N,-1,1), ncol=3, nrow=N)
##D y<-   10*exp( -rdist( x, rbind( c(.5,.5,.6) ) )/.5)
##D LKinfo<- LKrigSetup( x,  
##D                  LKGeometry = "LKBox",
##D                      nlevel = 1,
##D                      a.wght = 6.01,
##D                          NC = 20,
##D                   NC.buffer = 2,
##D                   normalize = FALSE,
##D              choleskyMemory = list(nnzR= 25e6),
##D                     lambda = .1,
##D                     mean.neighbor=200
##D                    ) 
##D # add in timing                    
##D  system.time(  out1<- LKrig( x,y, LKinfo=LKinfo) ) # about 25 seconds
##D # LKinfo$normalize<- TRUE
##D # system.time(  out2<- LatticeKrig( x,y, LKinfo=LKinfo) )# ~250 seconds
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
##D   # the Ring geometry will be periodic in the first dimension and rectagular on 
##D   # second. A useful approximation for spherical data omitting the polar caps. 
##D   
##D   LKinfo.CO2<- LKrigSetup(CO2$lon.lat, NC=100,nlevel=1, lambda=.2,
##D                        a.wght=5, alpha=1, 
##D                        LKGeometry="LKRing", choleskyMemory = list(nnzR=2e6) )
##D   print(LKinfo.CO2)                                          
##D   obj1<- LKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo.CO2)
##D # 5700+ basis functions 101X57 lattice  (number of basis functions
##D # reduced in y direction because of a rectangular domain
##D   obj1$trA.est # about 2900+ effective degrees of freedom 
##D #
##D   glist<- list( x= seq( -180,180,,200),y=seq( -80,80,,100) )
##D   fhat<- predictSurface( obj1,grid.list=glist)
##D #Plot data and gap-filled estimate
##D   set.panel(2,1)
##D   quilt.plot(CO2$lon.lat,CO2$y,zlim=c(373,381))
##D   title("Simulated CO2 satellite observations")
##D   world(add=TRUE,col="magenta")
##D   image.plot( fhat,zlim=c(373,381))
##D   world( add=TRUE, col="magenta")
##D   title("Gap-filled global predictions")
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
##D   LKinfo.test <- LKrigSetup( x, NC=16, nlevel=1, alpha=1.0,  a.wght=5)
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

## Not run: 
##D ##################################################################
##D #  a linear inverse problem 
##D ##################################################################
##D # NOTE the settings in the model are small to make for a fast example
##D #
##D data(ozone2)
##D x<- ozone2$lon.lat
##D y<- ozone2$y[16,]
##D good<- !is.na(y)	
##D x<- x[good,]
##D y<- y[good]
##D 
##D LKinfo<- LKrigSetup(x, a.wght=4.5, nlevel=3, nu=1, NC=4, lambda=.01)
##D 
##D # now observe a linear combination
##D   NNDist<- LKDist(x,x, delta=1.5) 
##D   A<- NNDist
##D   A$ra<- exp(-NNDist$ra)
##D # A is a weight matrix based on neighbors close by and
##D # in spind sparse matrix format
##D # now convert to spam format
##D   A<- spind2spam(A)
##D   
##D TMatrix <- get(LKinfo$fixedFunction)(x = x)
##D # Tmatrix is a 3 column matrix of constant and the two spatial coordinates
##D #  i.e. a linear function of the spatial variables 
##D PHI<- LKrig.basis(x, LKinfo)
##D X<-  A%*% PHI
##D U<-  A%*%TMatrix 
##D yIndirect<- A%*%y
##D #
##D # X<-  A##D 
##D 
##D out0<- LatticeKrig(x,y, LKinfo=LKinfo)
##D out1<- LatticeKrig(x,yIndirect, U=U, X=X, LKinfo=LKinfo)
##D 
##D # the predict function evaluates f in this case -- not the fitted values that
##D # are related to the 
##D # observations
##D # partial agreement but not exact -- this is because the observational models
##D # assume different errors
##D #
##D plot( predict(out0,x), predict( out1,x))
##D 
##D out<- LKrig.sim.conditional( out1,M=5, nx=10, ny=10 )
## End(Not run)	




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LKrig", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LKrig.basis")
### * LKrig.basis

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LKrig.basis
### Title: Functions for generating a multiresolution, compactly supported
###   basis, multiresolution covariance functions and simulating from these
###   processes.
### Aliases: LKrig.basis LKrig.precision LKrig.cov LKrig.cov.plot
###   LKrig.quadraticform LKrig.spind2spam LKrigCovWeightedObs
###   LKrigMarginalVariance
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
    nlevel <- 4
    a.wght <-  4 + 1/(.5)^2
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
##D     LKtemp<-  LKrigSetup( cbind( c( -1,1), c(-1,1)), NC=NC,
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
##D    LKtemp<-  LKrigSetup( cbind( c( -1,1), c(-1,1)), NC=NC,
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
##D   LKinfo<- LKrigSetup(cbind( c(-1,1), c(-1,1)), NC=6,
##D                                  nlevel=1, a.wght=4.5,alpha=1, NC.buffer=0 )
##D # evaluate the basis functions on a grid to look at them
##D   xg<- make.surface.grid( list(x=seq(-1,1,,50), y= seq(-1,1,,50)))
##D   PHI<- LKrig.basis( xg,LKinfo)
##D   dim(PHI) # should be  2500=50^2  by  36=6^2
##D # plot the 9th basis function  as.surface is a handy function to
##D # reformat the vector as an image object
##D # using the grid information in an attribute of the grid points
##D   set.panel(1,3) 
##D   image.plot(as.surface(xg, PHI[,9]))
##D   points(  make.surface.grid( LKrigLatticeCenters(LKinfo, 1)) , col="grey", cex=.5)
##D   title("A radial basis function")
##D # compare to the tensor product basis type
##D   LKinfo2<- LKrigSetup(cbind( c(-1,1), c(-1,1)), NC=6,
##D                                  nlevel=1, a.wght=4.5,alpha=1, NC.buffer=0,
##D                                  BasisType="Tensor" )
##D   PHI2<- LKrig.basis( xg,LKinfo2)
##D   image.plot(as.surface(xg, PHI2[,9]))
##D   points(  make.surface.grid( LKrigLatticeCenters(LKinfo, 1)), col="grey", cex=.5)
##D   title("Tensor product basis function")
##D   
##D   image.plot(as.surface(xg, PHI[,9] - PHI2[,9]))
##D   points(  make.surface.grid( LKrigLatticeCenters(LKinfo, 1)), col="grey", cex=.5)
##D   title(" Radial - Tensor for 9th basis function")                       
##D set.panel()
## End(Not run)
#
# example of basis function indexing
#
## Not run: 
##D # generating a basis on the domain [-1,1]X[-1,1] with 3 levels
##D # note that there are no buffering grid points.
##D   set.panel(3,2)
##D   LKinfo<-LKrigSetup(cbind( c(-1,1), c(-1,1)), NC=6,
##D                     a.wght=5, alpha=c(1,.5,.25), nlevel=3,
##D                     NC.buffer=0)
##D # evaluate the basis functions on a grid to look at them
##D   xtemp<- seq(-1,1,,40)
##D   xg<- make.surface.grid( list(x=xtemp, y= xtemp) )
##D   PHI<- LKrig.basis( xg,LKinfo)
##D # coerce to dense matrix format to make plotting easier.
##D   PHI<- spam2full(PHI)
##D # first tenth, and last basis function in each resolution level
##D # basis functions centers are added
##D  set.panel(3,3)
##D     for(  j in 1:3){
##D       id1<- LKinfo$latticeInfo$offset[j]+ 1
##D       id2<-  LKinfo$latticeInfo$offset[j]+ 10
##D       idlast<- LKinfo$latticeInfo$offset[j] +
##D                   LKinfo$latticeInfo$mx[j,1]*LKinfo$latticeInfo$mx[j,2]
##D    
##D       centers<-  make.surface.grid(LKrigLatticeCenters(LKinfo, j) )
##D       image.plot( as.surface(xg, PHI[,id1]))
##D       points( centers, cex=.2, col="grey")
##D       image.plot(as.surface(xg, PHI[,id2]))
##D       points( centers, cex=.2, col="grey")
##D       image.plot( as.surface(xg, PHI[,idlast]))
##D       points( centers, cex=.2, col="grey")
##D }
##D   set.panel()
## End(Not run)
## Not run: 
##D # examining the stationarity of covariance model
##D   temp.fun<- 
##D      function( NC.buffer=0, NC=4,  a.wght=4.01){
##D         LKinfo<- LKrigSetup(cbind( c(-1,1), c(-1,1)),nlevel=1, alpha=1,
##D                                  a.wght=a.wght, NC=NC,   
##D                                  NC.buffer=NC.buffer,
##D                                   choleskyMemory=list(nnzR=2e6))
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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LKrig.basis", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LKrig.sim")
### * LKrig.sim

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LKrig.sim
### Title: Functions for simulating a multiresolution process following the
###   Lattice Krig covariance model.
### Aliases: LKrig.sim LKrig.sim.conditional simConditionalDraw
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
  LKinfo<- LKrigSetup( x,NC=20,nlevel=1, alpha=1, lambda= .3 , a.wght=5)
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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LKrig.sim", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LKrigDefaultFixedFunction")
### * LKrigDefaultFixedFunction

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LKrigDefaultFixedFunction
### Title: Creates fixed part of spatial model.
### Aliases: LKrigDefaultFixedFunction LKrigPeriodicFixedFunction
###   predictLKrigFixedFunction
### Keywords: spatial

### ** Examples

x<- matrix( runif(100), nrow=50)
# linear polynomial 
T.matrix<- LKrigDefaultFixedFunction(x, m=2)
# quadratic polynomial 
T.matrix<- LKrigDefaultFixedFunction(x, m=3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LKrigDefaultFixedFunction", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LKrigLatticeCenters")
### * LKrigLatticeCenters

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LKrigLatticeCenters
### Title: Methods to report the locations or scales associated with the
###   lattice points.
### Aliases: LKrigLatticeCenters LKrigLatticeCenters.default
###   LKrigLatticeScales LKrigLatticeScales.default
###   LKrigLatticeCenters.LKBox LKrigLatticeCenters.LKInterval
###   LKrigLatticeCenters.LKRectangle LKrigLatticeCenters.LKRing
###   LKrigLatticeCenters.LKCylinder LKrigLatticeCenters.LKSphere
### Keywords: spatial

### ** Examples

  x<- cbind( c(-1,2), c(-1,2))
  LKinfo<- LKrigSetup( x, alpha=c( 1,.2,.01),
                   nlevel=3, a.wght=4.5, NC= 10)
# lattice centers for the second level   
# not points added for buffer outside of spatial domain                
   look<- LKrigLatticeCenters(LKinfo, Level=2) 
# convert grid format (gridList)  to just locations
   look<- make.surface.grid( look)
   plot( look,  cex=.5)
   rect( -1,-1,2,2, border="red4")
   
    x<- cbind( c(0, 360), c( 1,3)) 
    LKinfo<- LKrigSetup( x, LKGeometry="LKRing",
                   nlevel=1, a.wght=4.5, NC= 10, V= diag(c( 1,.01) ) )
                   
    polar2xy<- function(x){
	x[,2]*cbind( cos(pi*x[,1]/180), sin(pi*x[,1]/180))}
	        
    look1<- LKrigLatticeCenters( LKinfo, Level=1)               
    look2<- LKrigLatticeCenters( LKinfo, Level=1, physicalCoordinates=TRUE )
    look3<- polar2xy( look2$Locations )
# Basis scales:    
      LKrigLatticeScales( LKinfo)
    
    set.panel(3,1)
    plot( make.surface.grid( look1))
    plot( look2$Locations)
    plot( look3)

                 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LKrigLatticeCenters", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LKrigMiscellaneous")
### * LKrigMiscellaneous

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LKrig Miscellaneous Matrix Functions
### Title: Miscellaneous internal functions for LatticeKrig package.
### Aliases: LKrig.rowshift.periodic LKrig.shift.matrix LKArrayShift
###   LKrig.rowshift expandMList expandMatrix0 expandMatrix repMatrix
###   convertIndexPeriodic grid2Index convertIndexArray
### Keywords: spatial

### ** Examples

	A<- array( 1:90, c( 4,5,3))
	LKArrayShift( A, c( -1,-1,0))	
	
# welcome to the world of unrolling multiarray indices
A<- array( 1:60, c( 4,3,5))	
I<- rbind( c(1,2,1), c( 3,2,5))
look<- grid2Index( I, c( 4,3,5) )
# A has been filled with the right unrolled index
print( look)
print(A[look])



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LKrigMiscellaneous", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LKrigSAR")
### * LKrigSAR

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LKrigSAR
### Title: Method that creates the spatial autoregressive (SAR) matrix.
### Aliases: LKrigSAR LKrigSAR.LKBox LKrigSAR.LKRectangle
###   LKrigSAR.LKInterval LKrigSAR.LKRing LKrigSAR.LKCylinder
###   LKrigSAR.default LKrigSAR.LKSphere
### Keywords: spatial

### ** Examples

    x<- cbind( c(0,1))
	LKinfo<- LKrigSetup(x,LKGeometry="LKInterval",
	               nlevel=3, NC=3, a.wght=5, alpha=c(1,.5,.2) )
	B<- LKrigSAR( LKinfo, Level=2)
	B<-spind2full(B)
	image.plot( B)
	
	LKinfo<- LKrigSetup(cbind( c(0,360), c(0,1)) ,LKGeometry="LKRing",
	               nlevel=1, NC=3, a.wght=5, alpha=1)
	B<- LKrigSAR( LKinfo, Level=1)
	B<-spind2full(B)
	image.plot( B)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LKrigSAR", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LKrigSetup")
### * LKrigSetup

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LKrigSetup
### Title: Create or update the LatticeKrig model object (LKinfo) for
###   spatial fitting.
### Aliases: LKrigSetup LKinfoUpdate LatticeKrigEasyDefaults LKinfo
### Keywords: spatial

### ** Examples

  data(ozone2)
  # find the ranges of the  data, this is the same as passing
  # the entire set of observation locations and is more compact 
  sDomain<-apply( ozone2$lon.lat, 2,"range")
  LKinfo0<- LKrigSetup( sDomain, NC=10, nlevel=2, alpha=c(1,.5),
                       a.wght = 5)
  print( LKinfo0)
  
  #Gigantic buffer added note extra basis functions. 
  LKinfo<- LKrigSetup( sDomain, NC=10, NC.buffer= 15, nlevel=2, 
  alpha=c(1,.5),a.wght = 5)
  print( LKinfo)
  
  LKinfo2<- LKinfoUpdate( LKinfo,  a.wght=4.1, NC=12)
  LKinfo3<- LKrigSetup( sDomain, NC=12, nlevel=2, alpha=c(1,.5),
                        a.wght=4.1)
# LKinfo2 and LKinfo3 should be the same.                         
 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LKrigSetup", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LKrigSetupAlpha")
### * LKrigSetupAlpha

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LKrigSetupAlpha
### Title: Creates the alpha parameter list in LatticeKrig covariance.
### Aliases: LKrigSetupAlpha LKrigSetupAlpha.default
###   LKrigSetupAlpha.LKInterval LKrigSetupAlpha.LKRectangle
###   LKrigSetupAlpha.LKBox
### Keywords: spatial

### ** Examples

# an x that is just the limits of the domain	
  x<- cbind( c(0,1), c(0,1))
  
  LKinfo<- LKrigSetup( x, alpha=c( 1,.2,.01),
                   nlevel=3, a.wght=4.5, NC= 3)
  alphaList<- LKrigSetupAlpha( LKinfo)

  LKinfo<- LKrigSetup( x, nu=1, nlevel=4, a.wght=4.5, NC= 4)
  alphaList<- LKrigSetupAlpha( LKinfo)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LKrigSetupAlpha", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LKrigSetupAwght")
### * LKrigSetupAwght

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LKrigSetupAwght
### Title: Method to create a.wght component from the 'LKinfo' object.
### Aliases: LKrigSetupAwght LKrigSetupAwghtObject LKrigSetupAwght.default
###   LKrigSetupAwght.LKRectangle
### Keywords: spatial

### ** Examples

  x<- cbind( c(0,1))
  LKinfo<- LKrigSetup( x,LKGeometry="LKInterval", alpha=c( 1,.2,.01),
                   nlevel=3, a.wght=4.5, NC= 3)
  a.wghtList<- LKrigSetupAwght( LKinfo)
  
  x<- cbind( c(0,1), c(0,1))
  LKinfo<- LKrigSetup( x, alpha=c( 1,.2,.01),
                   nlevel=3, a.wght=4.5, NC= 3)
  a.wghtList<- LKrigSetupAwght( LKinfo)
# see   
  names(attributes( a.wghtList))
 
  



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LKrigSetupAwght", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LatticeKrig")
### * LatticeKrig

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LatticeKrig
### Title: User-friendly spatial prediction and inference using a compactly
###   supported multi-resolution basis and a lattice model for the basis
###   coefficients.
### Aliases: LatticeKrig print.LatticeKrig
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
##D # predict at observations:
##D   out.fhat<- predict( obj, xnew= x)
##D # coveniently predict on a 100X100 grid for plotting
##D  out.surf<- predictSurface( obj, nx=100, ny=100)
##D # image.plot( out.surf) 
## End(Not run)
# running an example by first setting up the model object
## Not run: 
##D # this is just a small model to run quickly
##D # compare the LKinfo object here  to one created implictily:  obj$LKinfo
##D LKinfo1<- LKrigSetup( x, NC=5, nlevel=3, a.wght=4.1, nu=1.0)
##D obj1<- LatticeKrig( x,y, LKinfo= LKinfo1)
## End(Not run)
#
# In this example lon/lat are treated as just Euclidean coordinates 
# a quick adjustment for small regions is to account for the difference
# in physical distance in N-S verses E_W
# is to just scale the longitude degrees to be comparable to degrees in latitude
# at least in the middle of the domain. The assumption is that for small spatial
# domains this approximation will not be bad for the coordinates at the edges too.
# You accomplish this by adding a scaling, V matrix:
# Here the V argument is rolled into the LKinfo object created within the function
#
## Not run: 
##D   meanLat<- mean( x[,2])*pi/180
##D   Vlonlat <- diag(  c( 1/cos(meanLat), 1) )
##D   obj1<- LatticeKrig( x, y, V = Vlonlat )
## End(Not run)

## Not run: 
##D # Refit using with just one level of  basis functions
##D # on a 20X20 grid witin the spatial domain ( so about 400) 
##D # actually number is 720 ( see obj1b$LKinfo) due adding edge nodes
##D # Add an aspect ratio of spatial domain 
##D # and find the a.wght parameter along with nugget and process variances.
##D # this takes a while partly because LatticeKrig model is not optimized for small data sets!
##D   obj1b<- LatticeKrig( x, y, nlevel=1, NC=20, findAwght=TRUE)
##D # rudimentary look at how likelihood was optimized
##D #log lambda and omega =  log(a.wght-4)/2 are useful parametrizations ...
##D   quilt.plot( obj1b$MLE$lnLike.eval[,c("logLambda","omega")],
##D        obj1b$MLE$lnLike.eval[,"lnProfileLike.FULL"], 
##D        xlab="loglamda", ylab="omega",
##D        zlim =c(-640,-612))
##D   points( obj1b$MLE$lnLike.eval[,c("logLambda","omega")],cex=.25)
##D       
## End(Not run)
# fitting replicate spatial data sets
# here we use the common observations over days for the ozone
# data set. Whether these are true replicated fields is in question
# but the analysis is still useful

## Not run: 
##D Y<-  na.omit( t( ozone2$y) ) 
##D ind<- attr( Y,"na.action")
##D X<- ozone2$lon.lat[-ind, ]
##D 
##D out1<- LatticeKrig( X, Y, nlevel=1, NC=20, findAwght=TRUE)
##D out2<- LatticeKrig( X, Y, nlevel=1, NC=20, findAwght=TRUE,
##D                         collapseFixedEffect=TRUE)
##D # compare the two models 
##D # Note second a.wght reflects more spatial correlation when individual 
##D # fixed effect is not removed ( 4.4 verses 4.07)
##D # nugget variance is nearly the same!
##D out1$MLE$summary[1:7]                        
##D out2$MLE$summary[1:7]
##D                         
##D 
##D 
## End(Not run)
## Not run: 
##D # Refit using the tensor product type of basis functions
##D # (default is "Radial"). An example how an additional argument that is 
##D # passed to the LKrigsetup function to create the LKinfo object.
##D   obj2<- LatticeKrig( x, y, BasisType="Tensor")
## End(Not run)

#
# A 1-d example with 3 levels of basis functions
# See LKrig for an explanation if nlevel, NC,  alpha and a.wght 
# covariance parameters.


## Not run: 
##D  x<- matrix(rat.diet$t)
##D  y<- rat.diet$trt
##D  fitObj<- LatticeKrig( x, y)
##D # NOTE lots of defaults are set for the model! See print( fitObj)
##D  plot( x,y)
##D  xg<- matrix(seq( 0,105,,100))
##D  lines( xg, predict(fitObj, xg) )
## End(Not run)

## Not run: 
##D #  a 3D example
##D set.seed( 123)
##D N<- 1000
##D x<-  matrix( runif(3* N,-1,1), ncol=3, nrow=N)
##D y<-   10*exp( -rdist( x, rbind( c(.5,.5,.6) ) )/.5)
##D 
##D # NOTE setting of memory size for Cholesky. This avoids some warnings and
##D # extra computation by the spam package
##D LKinfo<- LKrigSetup( x,  nlevel=1,  a.wght= 6.01, NC=6, NC.buffer=2,
##D                     LKGeometry="LKBox", normalize=FALSE, mean.neighbor=200,
##D                     choleskyMemory=list(nnzR= 2E6) )                                      
##D out1<- LatticeKrig( x,y, LKinfo=LKinfo)
##D 
##D glist<- list( x1=seq( -1,1,,30), x2=seq( -1,1,,30), x3 = 0)
##D xgrid<- make.surface.grid( glist)
##D 
##D yhat<- predict( out1, xgrid)
##D # compare yhat to true function created above
##D image.plot( as.surface( glist, yhat))
##D 
## End(Not run)
#
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
##D a.wghtGrid<-  4  +  c(.05, .1, .5, 1, 2, 4)^2
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
##D 
##D # MLE
##D obj0  <- LatticeKrig(CO.loc,  CO.tmin.MAM.climate, Z=CO.elev, 
##D                      findAwght=TRUE)
##D print(obj0$MLE$summary)
## End(Not run)

#########################################################################
# Reproducing some of the analysis for the example in the
# JCGS LatticeKrig paper.
#########################################################################

#### Here is an example of dealing with approximate spherical geometry.
## Not run: 
##D data(NorthAmericanRainfall)
##D library(mapproj)
##D x<- cbind(NorthAmericanRainfall$longitude, NorthAmericanRainfall$latitude)
##D y<- NorthAmericanRainfall$precip
##D log.y<- log(y)
##D elev<- NorthAmericanRainfall$elevation
##D # this is a simple projection as part of this and handled by the mapproj package
##D x.s<- mapproject( x[,1], x[,2], projection="stereographic")
##D x.s<- cbind( x.s$x, x.s$y)
##D 
##D # an alternative is to transform coordinates using another projection,
##D # e.g. a Lambert conformal projection
##D # with the project function from the rgdal package
##D # library( rgdal)
##D # x.s<- project(x,"+proj=lcc +lat_1=22 +lat_2=58 +lon_0=-93 +ellps=WGS84")
##D # this package has the advantage that the inverse projection is also 
##D # included ( inv=TRUE) so it is easy to evaluate the surface back on a Mercator grid.
##D              
##D obj0<- LatticeKrig(x.s, log.y, Z=elev )
##D 
##D fitSurface<- predictSurface( obj0, drop.Z=TRUE)
##D fitSurface$z<-  exp(fitSurface$z)/100
##D colorTable<- designer.colors( 256, c("red4", "orange", "yellow","green1", "green4", "blue"))
##D image.plot( fitSurface, col=colorTable)
##D map( "world", add=TRUE, col="grey30", lwd=3, proj="") 
##D 
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LatticeKrig", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PeriodicGeometry")
### * PeriodicGeometry

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Spatial models for data on spherical regions. 
### Title: Geometries to represent 2-d and 3-d spherical data.
### Aliases: LKRing LKCylinder LKSphere
### Keywords: spatial

### ** Examples

# 	
# fit the CO2 satellite data with a fixed lambda
# (but use a very small, unrealistic number of basis functions and levels so example
#  runs quickly)
 data(CO2)
# to look at raw data: quilt.plot(CO2$lon.lat, CO2$y, nx=288, ny=165)
# Use two levels starting at the second genertation of lattice points from the triangularization
  LKinfo0<- LKrigSetup( CO2$lon.lat, startingLevel=2 ,nlevel=2,
                       a.wght=1.1, alpha=c(1,.25),
                       LKGeometry="LKSphere" )
# Take a look at Model summary
  print( LKinfo0)
  
# Use arbitrary  lambda 
  obj0<- LKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo0, lambda=0.01)
  surface(obj0, nx=288, ny=165)
  world( add=TRUE)
## Not run: 
##D data(CO2)
##D # estimate lambda ( should be 0.003
##D # NOTE: lambda values will tend to be sensitive to the model choice
##D   LKinfo0<- LKrigSetup( CO2$lon.lat, startingLevel=2 ,nlevel=2,
##D                        a.wght=1.1, alpha=c(1,.25),
##D                        LKGeometry="LKSphere") 
##D   obj0B<-  LatticeKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo0)
##D # a more serious model this uses about 3300 basis functions
##D LKinfo0<- LKrigSetup( CO2$lon.lat, startingLevel=3, ,nlevel=3,
##D                        a.wght=1.1, alpha=c(1, .5, .25),
##D                        LKGeometry="LKSphere" )
##D                        
##D obj0B<-  LatticeKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo0)
##D # takes about 1 minute on a macbook air
##D # setting findAwght = TRUE  takes about 8 minutes with 
##D # lambda = 1.737 and a.wght = 15.8
## End(Not run) 
#####################################
# The ring geometry
#####################################
## Not run: 
##D   data(CO2)
##D   LKinfo1<- LKrigSetup(CO2$lon.lat, NC=8 ,nlevel=1, lambda=.2,
##D                        a.wght=5, alpha=1, 
##D                        LKGeometry="LKRing" )                                         
##D   obj1<- LatticeKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo1)	
##D # take a look: 
##D surface( obj1)
##D world( add=TRUE) 
## End(Not run)
# compare to fitting without wrapping:
## Not run: 
##D   LKinfo2<- LKrigSetup(CO2$lon.lat, NC=8 ,nlevel=1,
##D                    lambda=.2, a.wght=5, alpha=1 )                                         
##D   obj2<- LKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo2)	
##D  # NOTE: not periodic in longitude:
##D  surface(obj2)  
## End(Not run)

# a synthetic example and larger example
## Not run: 
##D  set.seed(124)
##D  N<- 1e4
##D   x0<- matrix( rnorm(3*N), ncol=3)
##D   x0<- x0/ sqrt( rowSums( x0^2))
##D   
##D   x<-  toSphere( x0 )
##D   
##D # the true function for testing -- a bump at the direction alpha
##D   fun<- function(X){
##D     alpha<-  c( .1,.1,1)
##D     alpha<- alpha/ sqrt( sum( alpha^2))
##D     4*( 1 + c(( X)%*%alpha) )^2 
##D   }
##D   
##D   ytrue <- fun(x0)
##D   y<- ytrue + .05*rnorm( length(ytrue))
##D # this defines about 3300 basis functions
##D   LKinfo1<- LKrigSetup( x,
##D                         startingLevel=3,
##D                         LKGeometry="LKSphere",
##D                         a.wght=1.01,
##D                         nlevel=3, alpha = c(1.0,.5,.25)^2,
##D                         choleskyMemory=list(nnzR= 20E6),
##D                         normalize=TRUE)
##D   out<- LatticeKrig( x,y, LKinfo=LKinfo1, lambda=.01)                      
##D surface( out)                        
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PeriodicGeometry", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Radial.Basis")
### * Radial.Basis

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Radial.basis
### Title: Two dimensional radial and tensor basis functions based on a
###   Wendland function.
### Aliases: Radial.basis LKrig.cyl WendlandFunction Tensor.basis triWeight
### Keywords: spatial

### ** Examples

set.seed(12)
x<- cbind( runif(100), runif(100))
center<- expand.grid( seq( 0,1,,5), seq(0,1,,5))
# coerce to matrix
center<- as.matrix(center)

  PHI1<- Radial.basis(x, center, basis.delta = .5)
  PHI2<- Tensor.basis( x, center, basis.delta = .5 )
# similarity of radial and tensor product forms  
  plot( c(0,1.1), c(0,1.1), type="p")
  for( k in 1:25){
	points( PHI1[,k], PHI2[,k])
	}
	
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
    
  obj1<- LKrig(x,y,NC=30,nlevel=1, alpha=1, 
    lambda=.01, a.wght=5, BasisFunction="triWeight", overlap=1.8)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Radial.Basis", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("directionCosines")
### * directionCosines

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: directionCosines
### Title: Utility functions for spherical coordinate and projections.
### Aliases: directionCosines toSphere projectionSphere
### Keywords: spatial

### ** Examples

# 
# icosohedron:
 phi = (1 + sqrt(5))/2
  V = matrix(c(0, 1, phi, 0, 1, -phi, 0, -1, phi, 0, -1, -phi, 
        1, phi, 0, -1, phi, 0, 1, -phi, 0, -1, -phi, 0, phi, 
        0, 1, -phi, 0, 1, phi, 0, -1, -phi, 0, -1), 12, 3, byrow = TRUE)

# check : library( rgl); plot3d( V, size=10, col="red4" )
# as lon lat:
V2<- toSphere( V)
plot( V2)

# lon lat grid
 lGrid<- make.surface.grid(  list(x= seq( -10,10,, 10), y= seq( -20,20,,10)) )
 
 dGrid<- directionCosines( lGrid)
 pairs( dGrid)
 # also try:   library( rgl); plot3d(  dGrid)
 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("directionCosines", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gridList-class")
### * gridList-class

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gridList-class
### Title: Class '"gridList"'. A description of a regular and
###   multidimensional grid.
### Aliases: gridList-class gridListInfo gridList
### Keywords: classes

### ** Examples

showClass("gridList")
# a 3-d grid
grid<- structure(
 list( x= seq( -1,1,,20), y= seq( 0,1,,15) ,oneMore = 1:10) ,
 class= "gridList" )



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gridList-class", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("nonstationaryModels")
### * nonstationaryModels

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: nonstationaryModels
### Title: Specifying nonstationary models
### Aliases: nonstationaryModels

### ** Examples


######################################################
##### This is an extended example showing how to define 
##### spatially varying rho parameter 
#######################################################
# Define some useful predict functions. 
#######################################################
  predict.surfaceGrid<- function(object,x){
    interp.surface( object, x)
    }
    
  predict.multivariateSurfaceGrid<- function(object,x){
    dimZ<- dim( object$z)
    L<- dimZ[3]
    out<- matrix( NA, nrow= nrow(x), ncol=L)
    for (  l in 1:L){
     out[,l]<- interp.surface( 
     list( x=object$x,y=object$y, z=object$z[,,l]) , x)
     }
     return( out)
  }
  
  predict.constantValue<- function(object,x){
   n<- length(object$values)
   m<- nrow( x)
   return( matrix( object$values, nrow=m, ncol=n, byrow=TRUE ) )
    }

################################################
##### Nonstationary examples
###############################################
# spatial domain    
sDomain<- rbind( c(-1.2,-1.2),
                 c(1,1))

# we will use this coarse grid to define any 
# surfaces of parameters
# (unrelated to the lattice grids and plotting grid!)
# this is larger than the sDomain to accomodate buffer points
# (with larger ranges when NC is small)
  gridList<- list( x = seq( -3, 3,,50),
                   y = seq( -3, 3,,75) )
  xPoints<- make.surface.grid( gridList)
  fineGrid<- make.surface.grid(
                 list( x = seq(-1, 1, ,45),
                       y = seq(-1, 1, ,60)
                       )
                               )
 
##################################################
### end of setup 
#################################################
# rho increases across the domain as a function of first coordinate. 
  rhoTemp<-  .01 +  10* pnorm( xPoints[,1], mean=.25, sd =.3 )
  rho.object<- as.surface( xPoints, rhoTemp) 
  class( rho.object)<- "surfaceGrid"
     
  LKinfo<- LKrigSetup( sDomain, NC= 4, nlevel = 3,
                 a.wght=4.5, nu=1, rho.object=rho.object)   
# simulate a field from this model
  set.seed(123)
  look<- LKrig.sim( fineGrid, LKinfo)
  image.plot( as.surface( fineGrid, look))
  xline( .25, col="grey30")
# see also 
# temp<- as.surface( fineGrid, look)
# I<- 20
# matplot(temp$x, temp$z, type="l", xlab="x", ylab="GP slice" )
  
######################################################
##### spatially varying alpha parameters 
#######################################################

# the alpha surface at each level will just be the result of 
# bilinear interpolation of values specified on a small grid.
# To keep things identified the alpha weights at 
# any grid location are 
# normalized to sum to 1. 
#
# create a 3 column matrix with  (proportional) alpha weights
# at each grid point 
#
  taper<- pnorm( xPoints[,1], mean = .4, sd=.02)
  alphaTemp<- cbind( taper,
        rep( 1, length( taper)), 
                        1-taper)
# normalize to sum to one                        
  alphaTemp <- alphaTemp/rowSums(alphaTemp)
 
# pack as a list 
# convert from a vector to the image/list format  $x $y $z
# give this object a class so that predict.surfaceGrid is called.
# accumulate these objects in a list 
# (yes this is a "list of lists")
  alphaObject<- list()
  for( k in 1:3){
     hold<- as.surface( xPoints, alphaTemp[,k]) 
     class( hold)<- "surfaceGrid"
     alphaObject<- c( alphaObject, list( hold))
  }
  
# define the 2-d LatticeKrig model
  LKinfo<- LKrigSetup(sDomain, NC = 4, a.wght=4.5,
              alpha = c(1,1,1), nlevel = 3, 
              alphaObject =  alphaObject )
# simulate a field 
 
  set.seed(123)
  look<- LKrig.sim( fineGrid, LKinfo)
  image.plot( as.surface( fineGrid, look))
  
######################################################
##### spatially varying a.wght parameters 
##### See above comments and setup
##### for steps that are the same 
#######################################################
  taper<- pnorm( xPoints[,1] + xPoints[,1],
                    mean = 0, sd=.15)
  a.wghtTemp<- 4.001*taper +  10*(1-taper)
# pack up as a list 
# convert from a vector to the image/list format  $x $y $z
# give this object a class so that predict.surfaceGrid is called.
# accumulate these objects in a list (yes this is 
# a "list of lists")
 
     a.wghtObjectA <- as.surface( xPoints, a.wghtTemp) 
     class( a.wghtObjectA)<- "surfaceGrid"
     

# define the 2-d LatticeKrig model
 
  LKinfo2<- LKrigSetup(sDomain, NC = 5, NC.buffer=0, 
              alpha = c(1, .5, .125), nlevel = 3, 
              a.wghtObject =  a.wghtObjectA)
              
  set.seed(123)            
  look<- LKrig.sim( fineGrid, LKinfo2)
  image.plot( as.surface( fineGrid, look))
##############################################
###### 1-d example
#############################################
  xCoarse1<- seq( -.5,1.5,, 40)
  y<-  pnorm( xCoarse1, mean=.4, sd=.05)*5 + 2.2 
  a.wghtObject<- Tps(xCoarse1, y, lambda=0)
  alphaTemp<-c(.5, .3, .2)
  LKinfoTEST<- LKrigSetup( rbind(0,1), NC=10,
                          LKGeometry="LKInterval",
                          nlevel=3, alpha=alphaTemp,
                          a.wghtObject = a.wghtObject,
                          NC.buffer=2
                          ) 
  xFine1<- cbind(seq( 0,1,length.out= 200))
  set.seed( 123)
  look<- LKrig.sim( xFine1, LKinfoTEST, M=5)
  matplot( xFine1, look, type="l", lty=1)
##################################################
######## Anisotropy in a.wght
##################################################
#### stationary example
  a.wghtM2<- c( rbind( c(  0,   0, -1.5),
                       c(-.5, 4.5,  -.25),
                       c(-1.5,  0,    0)
                   )
                   ) 
  
  LKinfo3<- LKrigSetup(sDomain, NC = 5, 
      a.wght= list( a.wghtM2), 
              alpha = c(1, .5, .125), nlevel = 3, 
              a.wghtObject =  NULL, normalize=TRUE )
              
  
  set.seed(123)            
  look<- LKrig.sim( fineGrid, LKinfo3)
  image.plot( as.surface( fineGrid, look))

#### Anisotropy varying over space
#### First check that the constant model can be reproduced
  a.wghtM2<- c( rbind( c(  0,   0, -1.5),
                       c(-.5, 4.5,  -.5),
                       c(-1.5,  0,    0)
                   )
                   )
                   
  a.wghtObject<- list( values=a.wghtM2)
  class(a.wghtObject )<- "constantValue"
 
  LKinfo4<- LKrigSetup(sDomain, NC = 5, 
              alpha = c(1,.5, .125), nlevel = 3, 
              a.wghtObject =  a.wghtObject, normalize=TRUE )
  set.seed(123)            
  look<- LKrig.sim( fineGrid, LKinfo4)
  image.plot( as.surface( fineGrid, look) )            

###### nonstationary anisotropy
 a.wghtA <- c( rbind( c(    0,   0, -2),
                       c( 0, 4.5,  0),
                       c(-2,  0,     0)
                   )
                   )
 a.wghtB <- c( rbind( c(  -2,   0,     0),
                       c(  0, 4.5,  0),
                       c(    0,    0, -2)
                   )
                   )
# Now create multivariate prediction object.
  gridList<-  attributes( xPoints)$grid.list
  m1<- length(gridList$x)
  m2<- length(gridList$y) 
  z<- array( NA, c( m1,m2,9))
  alpha<- (xPoints[,1] + 1 )/2
  alpha<- ifelse( alpha <= 0, 0, alpha)
  alpha<- ifelse( alpha >= 1, 1, alpha)
# A for loop over the 9 pixels   
  for(j in 1:9) {
# linear combination of two a.wght matrices
# 
    zTemp<- a.wghtA[j] * (1-alpha) +  a.wghtB[j]*(alpha)
# coerce into an image format    
    z[,,j]<- as.surface( xPoints, zTemp)$z
  }

  a.wghtObject<- list( x= gridList$x,  y= gridList$y, z=z )
  class( a.wghtObject)<- "multivariateSurfaceGrid"

  LKinfo5<- LKrigSetup(sDomain, NC = 25, NC.buffer=0,
              alpha = c(1,.5), 
              nlevel = 2, 
              a.wghtObject =  a.wghtObject )
  set.seed(122)  
  fineGrid<- make.surface.grid(
                 list( x = seq(-1, 1, ,150),
                       y = seq(-1, 1, ,180)
                       )
                               )
  look<- LKrig.sim( fineGrid, LKinfo5)
  image.plot( as.surface( fineGrid, look), col = terrain.colors(256) )




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("nonstationaryModels", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("registerdFORTRAN")
### * registerdFORTRAN

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: registeredFORTRAN
### Title: Internal FORTRAN routines for working with grids and finding
###   distances.
### Aliases: findnorm lkdist lkdistcomp lkdistgrid lkdistgridcomp
### Keywords: datasets

### ** Examples

print(lkdistgridcomp)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("registerdFORTRAN", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
