! this subroutine fills in the entries of the line/point distance matrix. The matrix is in a sparse format; entries are
! only populated when the closest distance between line/point pairs is less than the point's range.
! the number of entries must be known beforehand (via a call to LKTomGridCount) to make sure the right amount of space
! is allocated for the matrix entries.

! all the comments in this file (and any other file in this directory) that start with !$omp are OpenMP directives,
! which are instructions to the compiler that allow the program's work to be split into parallel threads.
subroutine LKTomGrid(dim, points, nPoints, lines, nLines, ranges, &
        rangeReps, nRanges, ind, entries, nEntries)
    implicit none
    integer, intent(in) :: dim, nPoints, nLines, nRanges, rangeReps(nRanges), nEntries
    integer :: rangeStarts(nRanges), outputIdx, localOutputIdx, rangeIdx, lineIdx, dimIdx, pointIdx
    integer, intent(out) :: ind(2, nEntries)

    double precision, intent(in) :: points(dim, nPoints), lines(2*dim, nLines), ranges(nRanges)
    double precision :: lineLengthSquared, dist, lineVec(dim), projectionResid(dim, dim), pointVec(dim), range
    double precision, intent(out) :: entries(nEntries)

    outputIdx = 1
    ind(:,:) = -1
    entries(:) = -5.1

    ! determining which points will have which ranges
    rangeStarts(1) = 0
    do rangeIdx = 1, (nRanges-1)
        rangeStarts(rangeIdx+1) = rangeStarts(rangeIdx) + rangeReps(rangeIdx)
    enddo
    !$omp parallel private(lineVec, lineLengthSquared, projectionResid, pointVec, dist, range)
        !$omp do
        do lineIdx = 1, nLines

            !making a matrix projectionResid from the line such that for a vector v, projectionResid times v is the
            !component of v perpendicular to the line
            lineVec = lines(:dim, lineIdx) - lines(dim+1:, lineIdx)
            lineLengthSquared = sum(lineVec * lineVec)
            projectionResid(:,:) = -1/lineLengthSquared
            do dimIdx = 1, dim
                projectionResid(dimIdx,:) = projectionResid(dimIdx,:) * lineVec(dimIdx)
                projectionResid(:,dimIdx) = projectionResid(:,dimIdx) * lineVec(dimIdx)
                projectionResid(dimIdx, dimIdx) = projectionResid(dimIdx, dimIdx) + 1
            enddo

            !iterate over the ranges and then the points with that range, so we don't need to get a new value for
            !the point's range every iteration
            do rangeIdx = 1, nRanges
                range = ranges(rangeIdx)
                do pointIdx = (rangeStarts(rangeIdx)+1), (rangeStarts(rangeIdx) + rangeReps(rangeIdx))
                    !using the matrix M from before to determine the distance between the line and point
                    !(divided by the range so the value is between 0 and 1)
                    pointVec = points(:, pointIdx) - lines(:dim, lineIdx)
                    pointVec = matmul(projectionResid, pointVec)
                    dist = sqrt(sum(pointVec * pointVec)) / range
                    if (dist < 1) then
                        !recording the distance, line number, and point number for the matrix entry
                        !$omp atomic capture
                            localOutputIdx = outputIdx
                            outputIdx = outputIdx + 1
                        !$omp end atomic
                        ind(1, localOutputIdx) = lineIdx
                        ind(2, localOutputIdx) = pointIdx
                        entries(localOutputIdx) = dist
                    endif
                enddo
            enddo
        enddo
        !$omp end do
    !$omp end parallel
end subroutine LKTomGrid