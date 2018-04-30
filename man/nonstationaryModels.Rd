\name{nonstationaryModels}
\alias{nonstationaryModels}
\title{Specifying a nonstationary model in alpha}
\description{An overview of specifying nonstationary models with some specific 2-d examples.}
\details{
The Lattice Krig model can be extended in a natural way to have a  nonstationary covariance that has a multiscale structure. 

The default and  process model has the form 

 \deqn{ g(x)= \sum_{l=1}^L g_l(x)
}
where g_l(x) are processes defined by fixed basis functions and with
coefficients that follow a mena zero multivariate normal distribution
and with dependence described by a spatial autoregrssion (SAR).

Under the model with approximate stationary the SAR is parametrized by
a set of a.wght values that are applied in the same way to every
lattice point and its neighbors. When the process is normalized to
have constant marginal variance, the variance of g_l(x) is given by
rho* alpha_l. We recommend that the alpha parameters sum to one and so
the marginal variance of g(x) is given by rho. This package allows for
the parameters a.wght, rho, and alpha to vary over the spatial
domain. In this way the variance of g_l(x) is given by
rho(x)alpha_l(x) and g(x) by rho(x) provided alpha_l sums to one. The
a.wght parameters can different values at each lattice point and
possibly at each level. Some preliminarly results suggest that having
a.wght vary differently for different levels may be too much
flexibility and may not needed.

We describe the formatting for these features in the LKinfo object below. Essentially they involve specifying one or more of the arguments: 
\code{a.wghtObject}, \code{rho.object} or \code{alphaObject}, in the call to 
\code{LKrigSetup}. 

\bold{a.wght}  The general form of this object is as a list of matrices. In this case
length( a.wght) = nlevel and a.wght[[l]] is a matrix where the number of rows is equal to the number of lattice points at level $l$ and in the order that they are indexed for the SAR matrix.  The number of columns depends on the particular geometry but we explain how this works for the rectangular spatial domain, LKRectangle. Here  a.wght[[l]]  has either 1 or 9 columns depending if anisotropy is specified. The weights for the LKRectangle anisotropy case are organized as
\preformatted{
    1 4 7
    2 5 8
    3 6 9
}
So an isotropic model with the center value as 4.5  looks like
\preformatted{
     0   -1    0
    -1   4.5  -1
     0   -1    0
}
and to encoded with the anisotropic extension the  row in 
a.wght[[l]] corresponding to this lattice point would be 

\code{ c( 0, -1, 0, -1, 4.5, -1, 0,-1, 0)}

Specifying the matrices of a.wghts directly can be invovled because
one needs to known the lattice information. An easier way to
accomplish specifying these models is to use a function on the spatial
domain to define the a.wght values based on the locations of the
lattice points. In this case one would pass an object, a.wghtObject,
that describes the function. This object needs to be defined so that
predict(a.wghtObject, x) will evaluate the function at the locations
x. (See example below.). With this object, LKrigSetupAwght will
evalute the function at the lattice poiints and so create the correct
list of matrices.

\bold{alpha}. To specify spatial varying alpha parameters one
specifies \code{alphaObject} as a list of objects where the predict
function works:

\code{predict(alphaObject[[l]], x) }

The result should give the values for alpha_l(x) with \code{x} a
matrix of arbitrary locations.

One should also set the vector of usual alpha parameters equal to all
 ones so only the spatially varying values as used for the variance.
 E.g. alpha= rep( 1, nlevel)
 

\bold{rho} The object \code{rho.object} is used define a spatially
varying rho(x). Like alpha and a.wght this object must work with the
predict function.

\code{predict(rho.object, x)}
}

\examples{

######################################################
##### This is an extended example showing how to define 
##### spatially varying rho parameter 
#######################################################
# Define a useful predict function. 
#######################################################
  predict.surfaceGrid<- function(object,x){
    interp.surface( object, x)
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
  a.wghtM2<- c( rbind( c(  0,   0, -1.5),
                       c(-.5, 4.5,  -.5),
                       c(-1.5,  0,    0)
                   )
                   ) 
  
  LKinfo3<- LKrigSetup(sDomain, NC = 5, 
      a.wght= list( a.wghtM2), 
              alpha = c(1,.5, .125), nlevel = 3, 
              a.wghtObject =  NULL, normalize=TRUE )
              
  
  set.seed(123)            
  look<- LKrig.sim( fineGrid, LKinfo3)
  image.plot( as.surface( fineGrid, look))
                  
}
  
