! this subroutine fills in the entries of the line/point distance matrix. The matrix is in a sparse format; entries are
! only populated when the closest distance between line/point pairs is less than the point's range.
! the number of entries must be known beforehand (via a call to LKTomGridCount) to make sure the right amount of space
! is allocated for the matrix entries.

! all the comments in this file (and any other file in this directory) that start with !$omp are OpenMP directives,
! which are instructions to the compiler that allow the program's work to be split into parallel threads.
subroutine LKTomGrid(dim, points, nPoints, lines, nLines, ranges, &
        rangeReps, nRanges, ind, entries, completions, nEntries)
    implicit none
    integer, intent(in) :: dim, nPoints, nLines, nRanges, rangeReps(nRanges), nEntries
    integer :: rangeStarts(nRanges), outputIdx, localOutputIdx, rangeIdx, lineIdx, dimIdx, pointIdx
    integer, intent(out) :: ind(2, nEntries)

    double precision, intent(in) :: points(dim, nPoints), lines(2*dim, nLines), ranges(nRanges)
    double precision :: dist2, lineVec(dim), pointVec(dim), range, range2, projDist, lineLength
    double precision, intent(out) :: entries(nEntries), completions(2,nEntries)

    outputIdx = 1
    ind(:,:) = -1
    entries(:) = -5.1
    completions(1,:) = 0
    completions(2,:) = 1

    ! determining which points will have which ranges
    rangeStarts(1) = 0
    do rangeIdx = 1, (nRanges-1)
        rangeStarts(rangeIdx+1) = rangeStarts(rangeIdx) + rangeReps(rangeIdx)
    enddo
    !$omp parallel private(lineVec, lineLength, pointVec, projDist, dist2, range, range2)
        !$omp do
        do lineIdx = 1, nLines

            !get the unit vector pointing along the line, which we will project the point onto
            lineVec = lines(dim+1:, lineIdx) - lines(:dim, lineIdx)
            lineLength = sqrt(sum(lineVec * lineVec))
            lineVec = lineVec / lineLength

            !iterate over the ranges and then the points with that range, so we don't need to get a new value for
            !the point's range every iteration
            do rangeIdx = 1, nRanges
                range = ranges(rangeIdx)
                range2 = range * range
                do pointIdx = (rangeStarts(rangeIdx)+1), (rangeStarts(rangeIdx) + rangeReps(rangeIdx))
                    !project the vector from the first line point to the point onto the line, then use the pythagorean
                    !theorem to get the squared distance from the point to the line
                    pointVec = points(:, pointIdx) - lines(:dim, lineIdx)
                    projDist = sum(lineVec * pointVec)
                    dist2 = sum(pointVec * pointVec) - projDist * projDist
                    if (dist2 < range2) then
                        !recording the distance, line number, and point number for the matrix entry
                        !$omp atomic capture
                            localOutputIdx = outputIdx
                            outputIdx = outputIdx + 1
                        !$omp end atomic
                        ind(1, localOutputIdx) = lineIdx
                        ind(2, localOutputIdx) = pointIdx
                        entries(localOutputIdx) = sqrt(dist2)/range
                            !calculating completions for each endpoint
                        if (sum(pointVec * pointVec) < range2) then
                            !the point is within range of the first endpoint
                                completions(1, localOutputIdx) = 0.5 - 0.5 * (projDist / range)
                        elseif (projDist < 0) then
                            !the point is before the start point, so the line integral is fully complete before it
                            completions(1, localOutputIdx) = 1
                        endif
                        projDist = projDist - lineLength
                        if (projDist**2 + dist2 < range2) then
                            !the point is within range of the second endpoint
                            completions(2, localOutputIdx) = 0.5 - 0.5 * (projDist / range)
                        elseif (projDist > 0) then
                            !the point is past the end point, so the line doesn't reach it
                            completions(2, localOutputIdx) = 0
                        endif
                    endif
                enddo
            enddo
        enddo
        !$omp end do
    !$omp end parallel
end subroutine LKTomGrid