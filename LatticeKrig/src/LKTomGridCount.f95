subroutine LKTomGridCount(dim, points, nPoints, lines, nLines, ranges, rangeReps, nRanges, nEntries)
    implicit none
    integer :: dim, nPoints, nLines, nRanges, rangeReps(nRanges)
    integer :: rangeStarts(nRanges), outputIdx, rangeIdx, lineIdx, dimIdx, pointIdx
    double precision :: points(dim, nPoints), lines(2*dim, nLines), ranges(nRanges)
    double precision :: lineLengthSquared, dist, lineVec(dim), projectionResid(dim, dim), pointVec(dim), range

    integer :: nEntries

    outputIdx = 1
    rangeStarts(1) = 0
    do rangeIdx = 1, (nRanges-1)
        rangeStarts(rangeIdx+1) = rangeStarts(rangeIdx) + rangeReps(rangeIdx)
    enddo
    outputIdx = 0
    !$omp parallel if(nLines > 50) private(lineVec, lineLengthSquared, projectionResid, pointVec, dist, range)
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