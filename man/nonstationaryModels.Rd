\name{nonstationaryModels}
\alias{nonstationaryModels}
\title{Specifying a nonstationary model in alpha}
\details{
The Lattice Krig model can be extended in a natural way to have a  nonstationary covariance that has a multiscale strucutre. 

The default and approximately stationary process model has the form 

 \deqn{ h(x)= \sum_{l=1}^L \sum_{j=1}^{m(j)} \phi_{m,l}(x) c_{m,l} } 
 
}
\examples{

######################################################
##### This is an extended example showing how to define 
##### spatially varying alpha parameters 
#######################################################

# the alpha surface at each level will just be the result of 
# bilinear interpolation of values specified on a small grid. 

# create a special purpose predict method to do this
# the tag surfaceGrid is the class name and of course can be more 
# specific for a project. 
  predict.surfaceGrid<- function(object,x){
    interp.surface( object, x)
   }
# spatial domain    
  sDomain<- rbind( c(-1,-1),
                     c(1,1))
# To keep things identified the alpha weights at 
# any grid location are 
# normalized to sum to 1. 
#
# coarse grid to define the alpha's 
# larger than spatial domain to accomdate for the buffer 
# points in lattice.
  
  gridList<- list( x = seq( -1.5,1.5,,20),
                   y = seq(-1.5,1.5,,25) )
  xPoints<- make.surface.grid( gridList)
  
# create a 3 column matrix with  (proportional) alpha weights
# at each grid point 
  taper<- pnorm( xPoints[,1], mean = .4, sd=.1)
  alphaTemp<- cbind( taper,
        rep( 1, length( taper)), 
                        1-taper)
  alphaTemp <- alphaTemp/rowSums(alphaTemp)
 
# pack as a list 
# convert from a vector to the image/list format  $x $y $z
# give this object a class so that predict.surfaceGrid is called.
# accumulate these objects in a list (yes this is a "list of lists")
  alphaObject<- list()
  for( k in 1:3){
     hold<- as.surface( xPoints, alphaTemp[,k]) 
     class( hold)<- "surfaceGrid"
     alphaObject<- c( alphaObject, list( hold))
  }
  
# define the 2-d LatticeKrig model
  LKInfo<- LKrigSetup(sDomain, NC = 4, a.wght=4.5,
              alpha = c(1,1,1), nlevel = 3, 
              alphaObject =  alphaObject )
# simulate a field 
  fineGrid<- make.surface.grid(
              list( x = seq( -1,1,,80),
                    y = seq( -1,1,,80))
                    )
  set.seed(123)
  look<- LKrig.sim( fineGrid, LKInfo)
  image.plot( as.surface( fineGrid, look))
  
######################################################
##### This is an extended example showing how to define 
##### spatially varying a.wght parameters 
##### See above comments for steps that are the same 
#######################################################
   sDomain<- rbind( c(-1,-1),
                     c(1,1))
# grid for evaluating a.wght has to be bigger in 
# order to include buffer points outside sDomain
   gridList<- list( x = seq( -1.8,1.8,,20),
                   y = seq(-1.8,1.8,,25) )
   xPoints<- make.surface.grid( gridList)
   
  taper<- pnorm( xPoints[,1] + xPoints[,1],
                    mean = 0, sd=.15)
  a.wghtTemp<- 4.001*taper +  10*(1-taper)
# pack up as a list 
# convert from a vector to the image/list format  $x $y $z
# give this object a class so that predict.surfaceGrid is called.
# accumulate these objects in a list (yes this is a "list of lists")
 
     a.wghtObject <- as.surface( xPoints, a.wghtTemp) 
     class( a.wghtObject)<- "surfaceGrid"
     
# define the 2-d LatticeKrig model
  LKInfo2<- LKrigSetup(sDomain, NC = 20, a.wght=4.5, 
              alpha = c(1,.5, .125), nlevel = 3, 
              a.wghtObject =  a.wghtObject,
              normalize = TRUE)
              
  fineGrid<- make.surface.grid(
              list( x = seq( -1,1,,100),
                    y = seq( -1,1,,100))
                    )
  set.seed(123)            
  look<- LKrig.sim( fineGrid, LKInfo2)
  image.plot( as.surface( fineGrid, look))
  
  
##################################################
######## Anisotropy in a.wght
##################################################
  a.wghtM2<- c( rbind( c(  0,   0, -1.5),
                       c(-.5, 4.5,  -.5),
                       c(-1.5,  0,    0)
                   )
                   ) 
  
  LKInfo3<- LKrigSetup(sDomain, NC = 20, 
      a.wght= list( a.wghtM2), 
              alpha = c(1,.5, .125), nlevel = 3, 
              a.wghtObject =  NULL, normalize=TRUE )
              
  fineGrid<- make.surface.grid(
              list( x = seq( -1,1,,100),
                    y = seq( -1,1,,100))
                    )
  set.seed(123)            
  look<- LKrig.sim( fineGrid, LKInfo3)
  image.plot( as.surface( fineGrid, look))
                  
}
  
