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


       subroutine dfind2d( x1, n1, x2, n2, delta2, ind, rd, Nmax, iflag)
       integer n1,n2,ind(Nmax,2)
       integer kk, i,j, ic
       real*8 x1(n1,2), x2(n2,2), delta2(n2), rd(Nmax), dtemp
c****   counter  for accumulating close points
        kk=0 
          do  15 i= 1, n1
            kksave= kk
            do 10 j =1,n2
c**** accumulate squared differences
              dtemp= (x1(i,1) - x2(j,1))**2 + (x1(i,2) - x2(j,2))**2
                if( dtemp.gt.delta2(j)) goto 10
c****       dtemp is less than D0 so save it as a close point
              kk=kk+1
c**** check if there is still array space 
              if( kk .gt. Nmax) then 
                iflag= -1
                return
              else
                ind(kk,1)= i
                ind(kk,2)= j
                rd(kk)= sqrt( dtemp)
              endif     
 10        continue          
 15      continue
         iflag=1
         Nmax=kk  
 20      continue
      end

  
      subroutine dfind3d( x1, n1, x2, n2, delta2, ind, rd, Nmax, iflag)
       integer n1,n2,ind(Nmax,2)
       integer kk, i,j, ic
       real*8 x1(n1,3), x2(n2,3), delta2(n2), rd(Nmax), dtemp
c****   counter  for accumulating close points
        kk=0 
          do  15 i= 1, n1
            do 10 j =1,n2
c**** accumulate squared differences
              dtemp= (x1(i,1) - x2(j,1))**2 + (x1(i,2) - x2(j,2))**2
     *               + (x1(i,3) - x2(j,3))**2   
                if( dtemp.gt.delta2(j)) goto 10
c****       dtemp is less than D0 so save it as a close point
              kk=kk+1
c**** check if there is still array space 
              if( kk .gt. Nmax) then 
                iflag= -1
                return
              else
                ind(kk,1)= i
                ind(kk,2)= j
                rd(kk)= sqrt( dtemp)
              endif     
 10        continue
 15      continue
         iflag=1
         Nmax=kk  
 20      continue
        end
 
