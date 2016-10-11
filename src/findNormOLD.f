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


C**   A fortran adaptation of findNorm()
c*** for the LatticeKrig package
c*** N. Lenssen October 2013 D Nychka November 15, 2013
          subroutine findNorm(mx,my,offset,Ux,Dx,Uy,Dy,
     *                        nLocations, xLocations, weights, Z)

	  integer mx, my, nLocations 
	  double precision offset
	  double precision Ux(mx,mx),Uy(my,my)
	  double precision Dx(mx),Dy(my)
	  double precision xLocations(nLocations,2)
          double precision weights(nLocations)
	  double precision Z(mx,my)
          integer iloc
          double precision normTemp
	  do iloc = 1, nLocations
             call findNormOne(mx,my,offset,Ux,Dx,Uy,Dy, 
     *                        xLocations(iloc,1), xLocations(iloc,2),
     *                        normTemp, Z )
             weights(iloc) = normTemp
          enddo
          return
          end

	  subroutine findNormOne(mx,my,offset,Ux,Dx,Uy,Dy,
     *                         kstar,lstar,normA,Z)
c* Note: Z is a work array.
	  integer mx,my
	  integer m1,n1,m2,n2
	  integer i,j,k,l
	  double precision offset, kStar, lStar, normA
	  double precision Ux(mx,mx),Uy(my,my)
	  double precision Dx(mx),Dy(my)	  
	  double precision Z(mx,my)
	  double precision distGrid, WendlandFunction
	  double precision normSum, ZUsum

	  normA = 0

c* find points in lattice where basis funtion evaluated at xi, xj may not
c* be zero.

	  m1 = max( ceiling(-offset+kStar),1)
	  n1 = max( ceiling(-offset+lStar),1)
	  m2 = min( floor(offset+kStar),mx)
	  n2 = min( floor(offset+lstar),my)	   

	  do k=m1, m2
	      do l=n1, n2
	          distGrid = sqrt( ((k-kStar)**2 + (l-lStar)**2))/offset
	      	  Z(k,l) = WendlandFunction(distGrid)
	      enddo
	  enddo

	  normSum = 0.0

	  do i=1, mx
	      do j=1, my
	   	      ZUsum = 0.0
	   	      do k=m1, m2
	   	          do l =n1, n2
	   	              ZUsum = ZUsum + Ux(k,i) * Uy(l,j) * Z(k,l)
	   	          enddo
	   	      enddo
	   	      normSum = normSum + (ZUsum / (Dx(i) + Dy(j)))**2
	   	  enddo
	  enddo

	  normA = normSum

	  return
	  end

	  double precision function WendlandFunction(d)
	  double precision d
	  if(d.LT.1) then
	      WendlandFunction = ((1 - d)**6*(35*d**2+18*d+3))/3
	  else
          WendlandFunction = 0
	  endif
	  return
	  end
 	
