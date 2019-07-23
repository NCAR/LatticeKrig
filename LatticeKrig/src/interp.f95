subroutine interp(grid, nGrid, delta, points, nPoints, output)
  integer nGrid, nPoints, gridIdx
  double precision delta, weight, grid(nGrid), points(nPoints)
  double precision output(nPoints)

  !$OMP PARALLEL NUM_THREADS(64) PRIVATE(gridIdx, weight)
    !$OMP DO
    do idx = 1, nPoints
      x = points(idx)
      gridIdx = 1 + FLOOR(x / delta)
      if (gridIdx >= nGrid) then
        output(idx) = grid(nGrid)
      else
        weight = x/delta - gridIdx + 1
        output(idx) = grid(gridIdx) * (1-weight) + grid(gridIdx + 1) * weight
      endif
    enddo
    !$OMP END DO
  !$OMP END PARALLEL

end
