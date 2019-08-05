!  this function does cubic spline interpolation using the provided grid, which is nGrid sets of 4 coefficients,
!  at points starting from 0 with separation delta. We compute the cubic spline interpolant at each of the nPoints
!  values in points, storing the results in output

subroutine CubicInterp(grid, nGrid, delta, points, nPoints, output)
    implicit none
    integer, intent(in) :: nGrid, nPoints
    integer :: gridIdx, idx

    double precision, intent(in) :: delta, grid(4, 0:nGrid-1), points(nPoints)
    double precision :: offset, x, coefs(4)
    double precision, intent(out) :: output(nPoints)

    !$OMP PARALLEL PRIVATE(gridIdx, offset, coefs)
        !$OMP DO SIMD
        do idx = 1, nPoints
            x = points(idx)
            gridIdx = FLOOR(x / delta)
            if (gridIdx >= nGrid) then
                ! the entry doesn't fit in the grid, so use the constant term of the last component of the cubic spline
                output(idx) = grid(1,nGrid-1)
            else
                ! grabbing the coeficients and evaluating the spline interpolant
                offset = x - gridIdx*delta
                coefs(:) = grid(:,gridIdx)
                output(idx) = coefs(1) + coefs(2)*offset + coefs(3)*offset*offset + coefs(4)*offset*offset*offset
            endif
        enddo
        !$OMP END DO SIMD
    !$OMP END PARALLEL
end subroutine CubicInterp