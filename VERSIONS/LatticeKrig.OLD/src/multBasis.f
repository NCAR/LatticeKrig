C** A fortran adaptation of findNorm()
c*** for the LatticeKrig package
c*** N. Lenssen October 2013 D Nychka November 15, 2013
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
	  integer i,j,k,l
	  double precision offset, kStar, lStar, normA
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

