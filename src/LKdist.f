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

       subroutine LKdist( x1, n1, x2, n2, dim,  delta2,
     *   ind, rd, Nmax, iflag)
       integer n1,n2, dim, Nmax, iflag
       integer ind(Nmax,dim)
       real*8 x1(n1,dim), x2(n2,dim), delta2, rd(Nmax)
       integer kk, i,j, l, ic
       real*8 dtemp
c****   counter  for accumulating close points
        kk=0 
          do  15 i= 1, n1
            kksave= kk
            do 10 j =1,n2
c**** accumulate squared differences
              dtemp =  0.0
              do 5 l = 1, dim
                dtemp= (x1(i,l) - x2(j,l))**2 + dtemp
 5            continue
                if( dtemp.gt.delta2 ) goto 10
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
      end
      
      subroutine LKdistComp( x1, n1, x2, n2, dim, delta,
     *  ind, rd, Nmax, iflag)
       integer n1,n2, dim, Nmax, iflag
       integer ind(Nmax,2)
       real*8 x1(n1,dim), x2(n2,dim), delta, rd(Nmax,dim)
       integer kk, i,j, l, ic
       real*8 dtemp(dim), dtest
c****   counter  for accumulating close points
        kk=0 
          do  15 i= 1, n1
            kksave= kk
            do 10 j =1,n2
c**** abs differences of each component
              do  l = 1, dim
                dtest = abs(x1(i,l) - x2(j,l))
                if( dtest.gt.delta ) goto 10
                dtemp(l)= dtest
              end do  
c****       all dtemp values are less than D0 so save it as a close point
              kk=kk+1
c**** check if there is still array space 
              if( kk .gt. Nmax) then 
                iflag= -1
                return
              else
                ind(kk,1)= i
                ind(kk,2)= j                
                do  l = 1, dim
                  rd(kk,l)=  dtemp(l)
                end do                  
              endif     
 10        continue          
 15      continue
         iflag=1
         Nmax=kk  
      end



 
