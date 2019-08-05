! this subroutine counts the number of entries the line/point distance matrix representation will need; that is, the
! number of line/point pairs such that the closest distance between them is less than the point's range.
! the calculations are all identical to the calculations in LKTomGrid.f95 (see that file for comments on the computation
! steps), the only differences are the parameter lists and what each subroutine does when it finds a line/point pair that
! is close enough; this file just increments outputIdx, while LKTomGrid records the distance, line number, and point
! number into the arrays that are passed back to R.

subroutine LKTomGridCount(dim, points, nPoints, lines, nLines, ranges, rangeReps, nRanges, nEntries)
    implicit none
    integer, intent(in) :: dim, nPoints, nLines, nRanges, rangeReps(nRanges)
    integer :: rangeStarts(nRanges), outputIdx, rangeIdx, lineIdx, dimIdx, pointIdx
    integer, intent(out) :: nEntries

    double precision, intent(in) :: points(dim, nPoints), lines(2*dim, nLines), ranges(nRanges)
    double precision :: lineLengthSquared, dist, lineVec(dim), projectionResid(dim, dim), pointVec(dim), range

    rangeStarts(1) = 0
    do rangeIdx = 1, (nRanges-1)
        rangeStarts(rangeIdx+1) = rangeStarts(rangeIdx) + rangeReps(rangeIdx)
    enddo

    outputIdx = 0
    !$omp parallel private(lineVec, lineLengthSquared, projectionResid, pointVec, dist, range)
    !$omp do reduction(+:outputIdx)
    do lineIdx = 1, nLines
        lineVec = lines(:dim, lineIdx) - lines(dim+1:, lineIdx)
        lineLengthSquared = sum(lineVec * lineVec)
        projectionResid(:,:) = -1/lineLengthSquared
        do dimIdx = 1, dim
            projectionResid(dimIdx,:) = projectionResid(dimIdx,:) * lineVec(dimIdx)
            projectionResid(:,dimIdx) = projectionResid(:,dimIdx) * lineVec(dimIdx)
            projectionResid(dimIdx, dimIdx) = projectionResid(dimIdx, dimIdx) + 1
        enddo

        do rangeIdx = 1, nRanges
            range = ranges(rangeIdx)*ranges(rangeIdx)
            do pointIdx = (rangeStarts(rangeIdx)+1), (rangeStarts(rangeIdx) + rangeReps(rangeIdx))
                pointVec = points(:, pointIdx) - lines(:dim, lineIdx)
                pointVec = matmul(projectionResid, pointVec)
                dist = sum(pointVec * pointVec)
                if (dist < range) then
                    outputIdx = outputIdx + 1
                endif
            enddo
        enddo
    enddo
    !$omp end do
    !$omp end parallel
    nEntries = outputIdx
end subroutine LKTomGridCount