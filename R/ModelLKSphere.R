
# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2012
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2

###################################################################
## LKrig model for data on Sphere using icosohedral grid.
##  Zach Thomas and Doug Nychka authors
####################################################################

###### Geometry class of this method is LKSphere

## These are obvious defaults for this model and saves
## specifying them in LKrigSetup
# This function is called in LKrigSetup before
# creating the lattice

setDefaultsLKinfo.LKSphere <- function(object, ...) {
# Definitely do not want Euclidean by default!
#  
  if( object$distance.type == "Euclidean"){
# Note default for lattice centers is on a unit sphere
# To change the radius in LKrigSetup pass the
# distance.type argument as "GreatCircle"
# but change its Radius attribute from 1.0.
      dType<- "GreatCircle"
      attr(dType, "Radius" ) <- 1.0
      object$distance.type <- dType
  }
# A lazy default: set alpha to 1.0 if only one level.
  if (object$nlevel == 1 & is.na(object$alpha[1])) {
    object$alpha <- list(1.0)
  }
#  
# hard wire the fixed part to just fit a constant function (m =1)
  if( !is.null( object$fixedFunction)){
    object$fixedFunction <- "LKrigDefaultFixedFunction"	
    object$fixedFunctionArgs$m <- 1
  }
# lazy default: set a.wght close to 6 giving a thin plate spline-like 
# model 
# (For the 12 points with 5 neighbors an adjustment is made
  if (is.na(object$a.wght)) {
    object$a.wght <- 6.01
  }
  return(object)
  }

# setup the lattice based on subdividing the faces of
# an icosohedron. There are 12 points at first
# level 
# Note that spatial coordiante passed in are assumed as 
# lon/lat in degrees with lon being [-180,180] to 
# be consistent with R maps package
 LKrigSetupLattice.LKSphere<- function(object, x=NULL, verbose,                                         
                                      NC=1, grid.info=NULL,
                                      ... ){
  if( is.null(x)){
    x<- object$x
  }
  ###### some common setup opertations to all geometries
  LKinfo<- object
  if(  class( LKinfo)[1] != "LKinfo") {
    stop("object needs to an LKinfo object")
  }
  
  rangeLocations<- apply( x,2, "range")
  nlevel<- LKinfo$nlevel
  ###### end common operations  
  # if x not passed then  assume full Sphere
  if( is.null(x)){
    x<- cbind( c(-90,90), c( -180,180))
  }
  if (is.null(NC)) {
    stop("Need to specify NC initial geodesic grid level")
  }
# Here NC is used as setting the initial or coarsest resolution geodesic grid. We may want to change the
# Right now, only allow use of
# up to 7th resolution grid due to memory issues
# if NC=3 and nlevel =5 then centers at resolutions 3,4, and 5 will
# be generated
  if (NC+nlevel-1 > 8){
    stop("NC+nlevel-1 cannot exceed 8")
  }else{
    ##This vector allows us to subset the full list of grids and just keep the ones we need.
    R<- seq(NC,NC+nlevel-1,1) 
    Rmax<- max(R)
  }
# delta cutoffs found empirically ...
# Nearest neighobors are within delta great circle distance (and second order
# neighbors are excluded)  
  delta<- 1.408/ 2^( 0:(Rmax-1) )
  delta.save<- delta[R]
  
##Build and subset geodesic grid up to level NC+nlevel-1; 
## returns each level in a list but in 
# 3-d coordinates (a.k.a. direction cosines)
  MultiGrid<- IcosohedronGrid(Rmax) ##Get full list of geodesic grids stopping at level Rmax
  grid.all.levels<- list()
  mLevel<- rep(NA,nlevel)
  for(l in (1:nlevel) ){
    # to Sphere converts to lon/lat coordinates
    gridTemp<- toSphere( MultiGrid[[ l + (NC -1) ]] )
      # trim to range of locations (in lon/lat coordinates)
    ind<-  gridTemp[,1] >= rangeLocations[1,1] &
           gridTemp[,1] <= rangeLocations[2,1] &
           gridTemp[,2] >= rangeLocations[1,2] &
           gridTemp[,2] <= rangeLocations[2,2] 
    ind2<-  (1:length( ind)) <= 12
# first 12 coordinates are always the initial isocosohedron points    
    grid.all.levels[[l]]<- gridTemp[ ind,] 
# logical attribute indicates which are initial vertices  
# within the subset determined by ranges     
    attr(grid.all.levels[[l]],"pentagon")<- ind2[ind]
    mLevel[l]<- nrow( grid.all.levels[[l]] )             
  }
  m<- sum(mLevel)
  offset <- as.integer(c(0, cumsum(mLevel)))
  out<- list( m=m,
              offset=offset,
              mLevel=mLevel, 
              delta=delta.save, 
              rangeLocations=rangeLocations,
# specific arguments for LKSphere              
              NC=NC,
              grid=grid.all.levels)
  return(out)
}

# LKrigSAR.LKSphere=function(
LKrigSAR.LKSphere = function(object, Level, ...) {
  if( Level>7){
    stop("can not handle more than 7 levels")}
# delta is the cutoff to find nearest neighbors -- tuned to this
# particular lattice.
  grid <- object$latticeInfo$grid[[Level]]
  delta <- object$latticeInfo$delta[[Level]]
  a.wght <- object$a.wght[[Level]]
## Find nearest 5 or 6 neighbors
  dType<- object$distance.type
  n <-  nrow(grid)
#  logical for pentagon points
# assume that the nodes are in same order as generated from
# IcosohedronGrid
  indP<- attr( grid,"pentagon")
  if( length( indP) != n) {
    print( c( n, length( indP) ))
    stop("Number of  elements in grid attribute different from grid")
  }
# sparse distance matrix in spind format.
# $ind are the indices that are nonzero
# Since original 12 iocosahedral vertices are stored as the first 12 points in each grid, we can easily
#  set their diagonal entries differently.
# The B matrix is already in spind format and the entries  (ra) are 
# converted to be the SAR matrix.
# The actual distances returned might be used to fine the weighting 
# as the distances vary some due the points not being on geodesics.
  B = LKDist(  grid[,1:2], grid[,1:2], delta = delta,
                  distance.type= dType)
# find the diagonal elements of the sparse matrix
  ind<- B$ind[,1] == B$ind[,2]
  if( sum( ind)!= nrow( grid)) {
    stop( "Number of diagonal elements in B different from grid")}
# take care of the 12 initial veritices with 5 NNs  
  a.wghtAdjusted<- ifelse( indP, a.wght, (a.wght - 6) + 5 )
  B$ra[ind ] <- a.wghtAdjusted
  B$ra[!ind] <- -1
  return(B)
# B is converted to spam format in LKrig.precision   
}

# return the lattice centers at a given Level
LKrigLatticeCenters.LKSphere<- function(object, Level, ... ){
  return( object$latticeInfo$grid[[Level]] )
} 

