module PolyIntFunctions
    private
    public :: MonomialIntegrals
    contains
    function MonomialIntegrals(L1, L2, maxDegree) result(integralValues)
        integer, intent(in) :: maxDegree
        integer :: pointIdx, entryIdx, m, n, mPlusN
        double precision, intent(in) :: L1(2), L2(2)
        double precision :: L(2), integralPoints(2,5)
        double precision :: integralValues((maxDegree+1)*(maxDegree+2)/2)
        double precision, parameter :: quadPoints(5) = (/ 0.046910077030668d0, 0.230765344947158d0, 0.500000000000000d0,&
                0.769234655052841d0, 0.953089922969332d0 /)
        double precision, parameter :: quadWeights(5) = (/ 0.118463442528095d0, 0.239314335249683d0, 0.284444444444444d0,&
                0.239314335249683d0, 0.118463442528095d0 /)

        L = L2 - L1
        do pointIdx = 1, 5
            integralPoints(:,pointIdx) = L1 + L*quadPoints(pointIdx)
        enddo
        entryIdx = 1
        do mPlusN = 0, maxDegree
            do n = 0, mPlusN
                m = mPlusN - n
                integralValues(entryIdx) = 1d0/(m+1) * sum(quadWeights * &
                        (integralPoints(1,:)**(m+1) * integralPoints(2,:)**n)) !x^(m+1) * y^n
                entryIdx = entryIdx + 1
            enddo
        enddo
        !L(2) is the x-component of the normal vector, pointing to the right from L1 to L2; we would use the unit vector,
        !but the numerical integration requires multiplying by the length of the line and getting the unit vector means
        !dividing by the length of the line, so these two cancel out
        integralValues = integralValues * L(2)
    end function MonomialIntegrals
end module PolyIntFunctions

subroutine PolyIntegral(vertices, nVerts, maxDegree, entries)
    use PolyIntFunctions
    integer, intent(in) :: nVerts, maxDegree
    integer :: vertIdx

    double precision, intent(in) :: vertices(2, nVerts)
    double precision :: L1(2), L2(2)
    double precision, intent(out) :: entries((maxDegree+1)*(maxDegree+2)/2)

    !evaluate the integral of monomials up to degree maxDegree using the divergence theorem; for each monomial x^m * y^n,
    !consider the vector-valued function <1/(m+1) x^(m+1) y^n, 0> which has divergence x^m * y^n, so integrating the
    !vector-valued function dotted with each edge's normal vector around the polygon will give us the integral of the
    !original monomial over the polygon
    entries(:) = 0
    !$omp parallel private(L1, L2)
    !$omp do reduction(+:entries)
    do vertIdx = 1, (nVerts-1)
        L1 = vertices(:,vertIdx)
        L2 = vertices(:,vertIdx+1)
        !if (vertIdx == nVerts) then
        !    L2 = vertices(:,1)
        !else
        !    L2 = vertices(:,vertIdx+1)
        !endif
        entries = entries + MonomialIntegrals(L1, L2, maxDegree)
    enddo
    !$omp end do
    !$omp end parallel
    entries = entries + MonomialIntegrals(vertices(:,nVerts), vertices(:,1), maxDegree)
end subroutine PolyIntegral