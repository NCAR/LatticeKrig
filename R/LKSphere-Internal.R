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
## Supporting function for LKSphere geometry for data on 
## Sphere using  an icosohedral grid.
##  Zach Thomas and Doug Nychka authors
####################################################################

# converts direction cosines to lon ( -180,180)  
# and lat coordinates in degrees 
toSphere=function(Grid){
  ##Convert to spherical in radians
  r=sqrt( rowSums(Grid^2) )
  lon=atan2(Grid[,2],Grid[,1])
  lat=asin(Grid[,3]/r)
  ##Output lon/lat in degrees
  LatLonGrid=cbind((180/pi)*lon,(180/pi)*lat)
  return(LatLonGrid)
  
}

projectionSphere<- function(x1, x2){
  x1<- rbind(x1)
  x2<- rbind(x2)
  # if on the north or south pole this is easy
  #
  if( abs(x1[1,3]) == 1.0 ){
    return( x2[,1:2] )
  }
  else{
    norm<- sqrt( x1[,1]^2 + x1[,2]^2)
    c1<- x1[,1]/norm
    s1<- x1[,2]/norm
    s2<- x1[,3]
    c2<- sqrt(1- s2^2)
    
    U<-  rbind( c( c1, s1, 0),
                c(-s1, c1, 0),
                c(  0,  0, 1))
    V<- rbind(  c( c2, 0, s2),
                c(  0, 1, 0),
                c(-s2,  0, c2))
    # the rotation matrix V%*%U rotates x1  around z to the prime meridan
    # 0 lon and then rotates in the x-z plane to the equator. 
    out<- V%*%(U%*%t(x2))
    
    return( t(out[2:3,]) )
  }
  
}


# generates multiresolution  grid to level K
IcosohedronGrid = function(K) {
  #require(abind)
  ##Golden Ratio
  phi = (1 + sqrt(5)) / 2
  ##Build Initial Icosahedron in cartesian coordinates
  V = matrix(
    c(
      0,   1, phi,
      0,   1,-phi,
      0,  -1, phi,
      0,  -1,-phi,
      1, phi,   0,
      -1, phi,   0,
      1,-phi,   0,
      -1,-phi,   0,
      phi,   0,   1,
      -phi,   0,   1,
      phi,   0,  -1,
      -phi,   0,  -1
    ),12,3,byrow = TRUE
  )
  V<- V/ sqrt( rowSums(V^2))
  MultiGrid = list(V)
#  
  if( K >1){  
# Define initial triangles using numerical reference labels and stored vertices
  Tri = array(rep(NA,20 * 3 * 3),dim = c(3,3,20))
# second place indexes the points, first the coordinates, third faces. 
  Tri[,,1] = t(V[c(1,3,9),])
  Tri[,,2] = t(V[c(1,3,10),])
  Tri[,,3] = t(V[c(3,8,10),])
  Tri[,,4] = t(V[c(3,7,9),])
  Tri[,,5] = t(V[c(1,5,9),])
  Tri[,,6] = t(V[c(1,6,10),])
  Tri[,,7] = t(V[c(1,5,6),])
  Tri[,,8] = t(V[c(3,7,8),])
  Tri[,,9] = t(V[c(6,10,12),])
  Tri[,,10] = t(V[c(8,10,12),])
  Tri[,,11] = t(V[c(7,9,11),])
  Tri[,,12] = t(V[c(5,9,11),])
  Tri[,,13] = t(V[c(2,5,11),])
  Tri[,,14] = t(V[c(2,5,6),])
  Tri[,,15] = t(V[c(2,6,12),])
  Tri[,,16] = t(V[c(4,8,12),])
  Tri[,,17] = t(V[c(4,7,8),])
  Tri[,,18] = t(V[c(4,7,11),])
  Tri[,,19] = t(V[c(2,4,12),])
  Tri[,,20] = t(V[c(2,4,11),])
#  for loop will add more levels to MuliGrid list
  for (i in 2:K) {
    Tri = aperm(Tri,c(2,1,3)) 
##Obtain bisection points for all three edges in each of the stored triangles
#
#                         A
#                      / (1) \
#                    AB   --  AC
#                  /   \ (3) /  \
#                B -(2)- BC - (4)- C
#    
    AB = (Tri[1,,] + Tri[2,,]) / 2
    BC = (Tri[2,,] + Tri[3,,]) / 2
    AC = (Tri[1,,] + Tri[3,,]) / 2
    ##There are many new vertices that are repeats...keep only the unique ones
    U = unique(t(cbind(AB,BC,AC)))
    ##Project onto sphere by normalizing vectors to radius one
    U <- U/sqrt(rowSums(U ^ 2))
    NewGrid = rbind(MultiGrid[[i - 1]],U)
    MultiGrid[[i]] = rbind(MultiGrid[[i - 1]],U)
    N = 20 * 4 ^ (i - 2) ##Number of triangles in new grid
    ##Make a list of all triangles for next increase in resolution. Four new trianges are formed within each
    ## of the old ones. Tri still has the old vertices...AB, AC, and BC have the midpoint vertices.
# Original code used abind, can avoid this because it is the last dimension being combined
#    Tri.n = array(rbind(Tri[1,,],AB,AC),dim = c(3,3,N))
#    Tri.n = abind(Tri.n, array(rbind(Tri[2,,], AB, BC),
#                               dim = c(3,3,N)), along = 3)
#    Tri.n = abind(Tri.n,  array(rbind(Tri[3,,],AC,BC),
#                               dim = c(3,3,N)),along = 3)
#    Tri.n = abind(Tri.n,  array(rbind(AB,BC,AC),
#                               dim = c(3,3,N)),along = 3)
    T1<- array(rbind(Tri[1,,],  AB,AC),
                dim = c(3,3,N))
    T2<- array(rbind(Tri[2,,], AB, BC),
                dim = c(3,3,N))
    T3<- array(rbind(Tri[3,,],AC,BC),
                dim = c(3,3,N))
    T4<- array(rbind(AB,BC,AC),
                dim = c(3,3,N))
    Tri.n<- array( c( T1,T2,T3,T4), c( 3,3, 4*N) ) 
#    test.for.zero( Tri.n, Tri.nTest)
    Tri = Tri.n
  }
}    
  return(MultiGrid)
}

# generates multiresolution  grid to level K
IcosohedronFaces = function(K) {
  #require(abind)
  ##Golden Ratio
  phi = (1 + sqrt(5)) / 2
  ##Build Initial Icosahedron in cartesian coordinates
  V = matrix(
    c(
      0,   1, phi,
      0,   1,-phi,
      0,  -1, phi,
      0,  -1,-phi,
      1, phi,   0,
      -1, phi,   0,
      1,-phi,   0,
      -1,-phi,   0,
      phi,   0,   1,
      -phi,   0,   1,
      phi,   0,  -1,
      -phi,   0,  -1
    ),12,3,byrow = TRUE
  )
  V<- V/ sqrt( rowSums(V^2))
  MultiGrid = list(V)
  #  
  if( K >1){  
    # Define initial triangles using numerical reference labels and stored vertices
    Tri = array(rep(NA,20 * 3 * 3),dim = c(3,3,20))
    # second place indexes the points, first the coordinates, third faces. 
    Tri[,,1] = t(V[c(1,3,9),])
    Tri[,,2] = t(V[c(1,3,10),])
    Tri[,,3] = t(V[c(3,8,10),])
    Tri[,,4] = t(V[c(3,7,9),])
    Tri[,,5] = t(V[c(1,5,9),])
    Tri[,,6] = t(V[c(1,6,10),])
    Tri[,,7] = t(V[c(1,5,6),])
    Tri[,,8] = t(V[c(3,7,8),])
    Tri[,,9] = t(V[c(6,10,12),])
    Tri[,,10] = t(V[c(8,10,12),])
    Tri[,,11] = t(V[c(7,9,11),])
    Tri[,,12] = t(V[c(5,9,11),])
    Tri[,,13] = t(V[c(2,5,11),])
    Tri[,,14] = t(V[c(2,5,6),])
    Tri[,,15] = t(V[c(2,6,12),])
    Tri[,,16] = t(V[c(4,8,12),])
    Tri[,,17] = t(V[c(4,7,8),])
    Tri[,,18] = t(V[c(4,7,11),])
    Tri[,,19] = t(V[c(2,4,12),])
    Tri[,,20] = t(V[c(2,4,11),])
    #  for loop will add more levels to MuliGrid list
    TriList<-list( 1:(K-1) )
    for (i in 2:K) {
      Tri = aperm(Tri,c(2,1,3)) 
      TriList[[i-1]]<- Tri
      ##Obtain bisection points for all three edges in each of the stored triangles
      #
      #                         A
      #                      / (1) \
      #                    AB   --  AC
      #                  /   \ (3) /  \
      #                B -(2)- BC - (4)- C
      #    
      AB = (Tri[1,,] + Tri[2,,]) / 2
      BC = (Tri[2,,] + Tri[3,,]) / 2
      AC = (Tri[1,,] + Tri[3,,]) / 2
      ##There are many new vertices that are repeats...keep only the unique ones
      U = unique(t(cbind(AB,BC,AC)))
      ##Project onto sphere by normalizing vectors to radius one
      U <- U/sqrt(rowSums(U ^ 2))
      NewGrid = rbind(MultiGrid[[i - 1]],U)
      MultiGrid[[i]] = rbind(MultiGrid[[i - 1]],U)
     
      N = 20 * 4 ^ (i - 2) ##Number of triangles in new grid
      ##Make a list of all triangles for next increase in resolution. Four new trianges are formed within each
      ## of the old ones. Tri still has the old vertices...AB, AC, and BC have the midpoint vertices.
      
      T1<- array(rbind(Tri[1,,],  AB,AC),
                 dim = c(3,3,N))
      T2<- array(rbind(Tri[2,,], AB, BC),
                 dim = c(3,3,N))
      T3<- array(rbind(Tri[3,,],AC,BC),
                 dim = c(3,3,N))
      T4<- array(rbind(AB,BC,AC),
                 dim = c(3,3,N))
      Tri.n<- array( c( T1,T2,T3,T4), c( 3,3, 4*N) ) 
      #    test.for.zero( Tri.n, Tri.nTest)
      # project triangles onto unit sphere note here columns index the vertices
      for ( k in 1: dim( Tri.n)[3]){
        for( l in 1:3){
        Tri.n[,l,k]<-Tri.n[,l,k]/ sqrt( sum( Tri.n[,l,k]^2 ) )
        }
      }
      
      Tri = Tri.n
    }
  }    
  # normlist points defining faces
  
  return(list(nodes=MultiGrid, Faces=TriList))
}
