subroutine CubicInterp(grid, nGrid, delta, points, nPoints, output)
    integer, intent(in) :: nGrid, nPoints
    double precision, intent(in) :: delta, grid(4, 0:nGrid-1), points(nPoints)

    integer :: gridIdx
    double precision :: offset, x, coefs(4)

    double precision, intent(out) :: output(nPoints)
    !$OMP PARALLEL PRIVATE(gridIdx, offset, coefs)
        !$OMP DO SIMD
        do idx = 1, nPoints
            x = points(idx)
            gridIdx = FLOOR(x / delta)
            if (gridIdx >= nGrid) then
                output(idx) = 0
            else
                offset = x - gridIdx*delta
                coefs(:) = grid(:,gridIdx)
                output(idx) = coefs(1) + coefs(2)*offset + coefs(3)*offset*offset + coefs(4)*offset*offset*offset
            endif
        enddo
        !$OMP END DO SIMD
    !$OMP END PARALLEL
end subroutine CubicInterp