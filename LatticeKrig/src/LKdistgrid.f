c****  # LatticeKrig  is a package for analysis of spatial data written for
c****  # the R software environment .
c****  # Copyright (C) 2016
c****  # University Corporation for Atmospheric Research (UCAR)
c****  # Contact: Douglas Nychka, nychka@ucar.edu,
c****  # National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
c****  #
c****  # This program is free software; you can redistribute it and/or modify
c****  # it under the terms of the GNU General Public License as published by
c****  # the Free Software Foundation; either version 2 of the License, or
c****  # (at your option) any later version.
c****  # This program is distributed in the hope that it will be useful,
c****  # but WITHOUT ANY WARRANTY; without even the implied warranty of
c****  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c****  # GNU General Public License for more details.


c     234567
c     Finds all distances in x1 that are within delta of a set of grid points.
c     results are returned in the row, column,  value sparse matrix format
c      
c     x1:  coordinates  n1 by ndim with columns being dimension and row indexing the
c     cases
c    
c     nGrid: gives the sizes of the grid in each dimension
c     and is it assumed each grid starts at 1 and the points are integer spaced.
c     delta cutoff distance only distances within delta are found
c     matrix of row and column indices
c     RETURNED ARGUMENTS      
c     irow: row indices of points within delta  (from 1 to n1)
c     icol: col indices of points wihtin delta   (from 1 to total number of grid points -- product of nGrid).
c     ra: distance
c     Nmax: maximum number of points allocated for storage. On return number of elements found that are non zero.
c     iflag: if  -1 then an error.
c INDEXING GRID:
c rightmost dimension index varies fastest
c for a  grid of  x= 1:5, y=1:4  z = 1:3 there are 5*4*3 = 60 total points
c the points are ordered as   expand.grid( 1:5, 1:4, 1:3) = 
c   1,1,1
c   2,1,1 
c    ...
c   5,1,1 
c   1,2,1
c   2,2,1
c    ...
c   4,4,3
c   5,4,3 
c    NOTE: This routine assumes that the coordinates have been scaled to 
c           to a unit grid spacing. So delta also is in the units where 
c           grid points have a unit spacing. 
c     
      subroutine lkdistgrid( x1, n1, nGrid, nDim, delta,
     *     irow, jcol, ra, Nmax, iflag)
      integer n1, nDim, Nmax
      double precision x1(n1,nDim)
      integer nGrid(nDim)
      double precision delta
      integer irow(Nmax), jcol(Nmax)
      double precision ra(Nmax)
      integer iflag
      double precision delta2, deltaX(10)
      double precision dist2Temp, x1ijTemp
      integer gridStep(10), cM(10), M(10)
      integer L
      integer indTemp, k, kk, kI, offset, prodM
      integer m1, m2
c
      delta2 = ( delta)**2
c work arrays only dimensioned up to 10 dimesions -- unlikely to be more 
c than 4       
      if( nDim.gt.10) then 
          iflag = -1
          Nmax = 0
          return
      endif    
c gridStep is the stride used to index each dimension
      gridStep(1) = 1
      do  j =  2, nDim
        gridStep(j) =  gridStep(j-1)* nGrid(j-1)
      enddo
c begin loop over coordinate points (x1) 
c kk is the counter for grid points that are within delta of the x1's   
      kk =  0
      do  i = 1, n1 
        prodM  =  1
        offset  =  0 
        do  j = 1, nDim 
          x1ijTemp = x1(i,j)
          m1   =  max( ceiling( -delta + x1ijTemp), 1 )
          m2   =  min( floor( delta + x1ijTemp), nGrid(j)) 
c jump to next case if this point has a coordinate that is delta away from box edge.           
          if( (m1.gt.nGrid(j)) .or. ( m2.lt.1) ) then
            goto 100
          endif         
c jth dimension has grid points from m1:m2 within delta of  
c   x(i,j). This is with the edges taken into account. 
c there are M(j) points in this set  
          M(j)  =  ( m2 - m1 ) + 1          
c prodM is the total number of grid points within a delta box of coordinate
c    x(i,j).          
          prodM =  prodM * M(j)
c    tranlate coordinates of x1(i,)  relative to the lower corner of delta box
          deltaX(j) =  x1ijTemp - m1
c offset is the position in the grid array to index the lower corner of the delta
c box  i.e. this corner point is offset + 1
          offset  =  offset + ( m1 - 1  )*gridStep(j)         
        enddo
        cM(1) =  1
        do  j = 2, nDim 
          cM(j) =  cM(j-1) *  M(j-1) 
        enddo
c     loop over all  grid points in delta box
c     unusual indexing is to allow a single to
c     handle a multidimensional grid.
        do  k = 1, prodM
c for a given grid point (k) loop over dimensions and accumulate
c the distance between it and the value for the point.
c also accumulate the index of this grid point for the 
c sparse distance matrix.
c kI is the grid point coordinate at dimension j
c minus its index at the "lower corner" of the delta box           
          indTemp =  0 
          dist2Temp =  0.0D0
          L =  k - 1   
          do  j = nDim, 1, -1 
        	kI =   L/cM(j)  
            indTemp =  indTemp +  kI * gridStep(j)
            dist2Temp =  dist2Temp + (deltaX(j) - kI)**2
            L  =  L - kI*cM(j)
          enddo 
c save distance if radial distance is <= delta            	
          if( dist2Temp .le. delta2) then
             kk =  kk + 1
c does kk overrun the size of the allocated space?             
             if( kk. gt. Nmax) then
                iflag = -1
                return
             endif
c the Euclidean distance              
             ra(kk)   =  dsqrt( dist2Temp)
c i^th point             
             irow(kk) =  i
c index of the grid point             
             jcol(kk) =  indTemp + 1 + offset
          endif
        enddo 
c next point         
 100  continue         
      enddo
c on return Nmax is actual number of entries in distance matrix      
      Nmax =  kk
c all is well ...      
      iflag=0
      end

      subroutine lkdistgridcomp( x1, n1, nGrid, nDim, delta,
     *     irow, jcol, ra, Nmax, iflag)
      integer n1, nDim, Nmax
      double precision x1(n1,nDim)
      integer nGrid(nDim)
      double precision delta
      integer irow(Nmax), jcol(Nmax)
      double precision ra(Nmax,nDim)
      integer iflag
      double precision delta2, deltaX(10), x1ijTemp
      double precision distComp(10), distTemp
      integer gridStep(10), cM(10), M(10)
      integer L
      integer indTemp, k, kk, kI, offset, prodM
      integer m1, m2
c
      delta2 = ( delta)**2
c work arrays only dimensioned up to 10 dimesions -- unlikely to be more 
c than 4       
      if( nDim.gt.10) then 
          iflag = -1
          Nmax = 0
          return
      endif    
c gridStep is the stride used to index each dimension
      gridStep(1) = 1
      do  j =  2, nDim
        gridStep(j) =  gridStep(j-1)* nGrid(j-1)
      enddo
c begin loop over coordinate points (x1) 
c kk is the counter for grid points that are within delta of the x1's   
      kk =  0
      do  i = 1, n1 
        prodM  =  1
        offset  =  0 
        do  j = 1, nDim 
          x1ijTemp = x1(i,j)
          m1   =  max( ceiling( -delta + x1ijTemp), 1 )
          m2   =  min( floor( delta + x1ijTemp), nGrid(j)) 
c jump to next case if this point has a coordinate that is delta away from box edge.           
          if( (m1.gt.nGrid(j)) .or. ( m2.lt.1) ) then
            goto 100
          endif         
c jth dimension has grid points from m1:m2 within delta of  
c   x(i,j). This is with the edges taken into account. 
c there are M(j) points in this set  
          M(j)  =  ( m2 - m1 ) + 1          
c prodM is the total number of grid points within a delta box of coordinate
c    x(i,j).          
          prodM =  prodM*(M(j))
c    tranlate coordinates of x1(i,)  relative to the lower corner of delta box
          deltaX(j) =  x1ijTemp - m1
c offset is the position in the grid array to index the lower corner of the delta
c box  i.e. this corner point is offset + 1
          offset  =  offset + ( m1 - 1  )*gridStep(j)         
        enddo
        cM(1) =  1
        do  j = 2, nDim 
          cM(j) =  cM(j-1)*( M(j-1))
        enddo
c loop over all points in delta box        
        do  k = 1, prodM
c for a given grid point loop over dimensions and accumulate
c the distance between it and the value for the point.
c also accumulate the index of this grid point for the 
c sparse distance matrix.
c kI is the grid point coordinate at dimension j
c minus the value at the "lower corner" of the delta box           
          indTemp =  0 
          L =  k - 1   
          do  j = nDim, 1, -1 
        	kI =   L/cM(j)  
            indTemp =  indTemp +  kI * gridStep(j)
            distTemp =  abs( deltaX(j) - kI)
            L  =  L - kI*cM(j)
            if(distTemp .ge. delta ) then 
              goto 90   
            endif
            distComp(j)= distTemp
          enddo 
c save distance if all component distances are <= delta            	
          kk =  kk + 1
c does kk overrun the size of the allocated space?             
          if( kk. gt. Nmax) then
                iflag = -1
                return
          endif
c copy distance components to output array
          do j = 1, nDim             
             ra(kk,j)   =  distComp(j)
          enddo
c i^th point             
             irow(kk) =  i
c index of the grid point             
             jcol(kk) =  indTemp + 1 + offset
c GOTO point if component is greater than delta  
 90     continue        
        enddo        
c GOTO to skip this point and try next one
 100  continue         
      enddo
c on return Nmax is actual number of entries in distance matrix      
      Nmax =  kk
c all is well ...      
      iflag=0
      end
      
c234567
      subroutine lkdistgrid2( x1, nx1, mx, my, delta2, ind,
     + rd, Nmax, iflag)
      double precision x1(nx1,2)
      integer nx1
      integer mx,my 
      double precision delta2
      integer ind(Nmax,2)
      double precision  rd(Nmax)
      integer Nmax
      integer iflag
      double precision kstar, lstar
      integer m1,n1,m2,n2
      integer kk, i
      double precision  delta, dtemp
c****   counter   do  accumulating close points
c234567
      delta = sqrt( delta2)
      kk = 0 
      do i = 1, nx1
        kStar = x1(i,1)      
        lStar = x1(i,2)            
        m1 = max( ceiling(-delta+kStar), 1)
        n1 = max( ceiling(-delta+lStar), 1)
        m2 = min( floor(   delta+kStar),mx)
        n2 = min( floor(   delta+lStar),my)
c**** loop over all grid points that are in +- delta square about (kstar, lstar)
c**** excluding edges 	  	   
        do k=m1, m2
           do l=n1, n2
             dtemp =  ((k-kStar)**2 + (l-lStar)**2)
             if( dtemp .le. delta2) then
               kk=kk+1
c****         check if there is still array space 
               if( kk .gt. Nmax) then 
                    iflag= -1
                    return
                 else
                    ind(kk,1)= i
                    ind(kk,2)= k + (l-1)*my
                    rd(kk)= sqrt(dtemp)
                 endif 
              endif
           enddo
        enddo
c23456
      enddo 
      Nmax=kk 
      iflag=1 
      end

 
