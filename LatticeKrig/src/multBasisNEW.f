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

          subroutine multBasis(mx,my,offset,CoefLevel,
     *                              nLocations, xLocations,f)
          integer mx, my, nLocations 
          double precision offset
          double precision CoefLevel(mx,my)
          double precision xLocations(nLocations,2)
          double precision floc, f(nLocations)
          integer iloc
          do iloc = 1, nLocations
             call multBasisOne(mx, my, offset, CoefLevel, 
     *                      xLocations(iloc,1), xLocations(iloc,2),floc)
             f(iloc) = floc
          enddo
          return
          end

          subroutine  multBasisOne(mx, my, offset, CoefLevel,
     *                                      kStar, lStar,floc)
          integer k,l
          double precision offset, kStar, lStar
          double precision CoefLevel(mx,my)
          double precision distGrid, WendlandFunction
          double precision floc
          double precision sumTerms
c* find points in lattice where basis funtion evaluated at xi, xj may not
c* be zero.
          m1 = max( ceiling(-offset+kStar),1)
          n1 = max( ceiling(-offset+lStar),1)
          m2 = min( ceiling( offset+kStar),mx)
          n2 = min( ceiling( offset+lstar),my)     
          sumTerms = 0
          do k=m1, m2
              do l=n1, n2
                distGrid = sqrt( ((k-kStar)**2 + (l-lStar)**2))/offset
                sumTerms =  sumTerms
     *                  + WendlandFunction(distGrid)* CoefLevel(k,l)
              enddo
          enddo
          floc = sumTerms
          return
          end

