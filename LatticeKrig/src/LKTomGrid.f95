subroutine LKTomGrid(dim, points, nPoints, lines, nLines, ranges, &
        rangeReps, nRanges, ind, entries, nEntries)
    implicit none
    integer :: dim, nPoints, nLines, nRanges, rangeReps(nRanges)
    integer :: nEntries
    integer :: rangeStarts(nRanges), outputIdx, rangeIdx, lineIdx, dimIdx, pointIdx
    double precision :: points(dim, nPoints), lines(2*dim, nLines), ranges(nRanges)
    double precision :: lineLengthSquared, dist, lineVec(dim), projectionResid(dim, dim), pointVec(dim), range

    integer :: ind(2, nEntries)
    double precision :: entries(nEntries)

    outputIdx = 1
    ind(:,:) = -1
    entries(:) = -5.1
    rangeStarts(1) = 0
    do rangeIdx = 1, (nRanges-1)
        rangeStarts(rangeIdx+1) = rangeStarts(rangeIdx) + rangeReps(rangeIdx)
    enddo
    !$omp parallel if(nLines > 50) private(lineVec, lineLengthSquared, projectionResid, pointVec, dist, range)
        !$omp do
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
                range = ranges(rangeIdx)
                do pointIdx = (rangeStarts(rangeIdx)+1), (rangeStarts(rangeIdx) + rangeReps(rangeIdx))
                    pointVec = points(:, pointIdx) - lines(:dim, lineIdx)
                    pointVec = matmul(projectionResid, pointVec)
                    dist = sqrt(sum(pointVec * pointVec)) / range
                    if (dist < 1) then
                        !$omp critical(critical_writeOutput)
                            ind(1, outputIdx) = lineIdx
                            ind(2, outputIdx) = pointIdx
                            entries(outputIdx) = dist
                            outputIdx = outputIdx + 1
                        !$omp end critical(critical_writeOutput)
                    endif
                enddo
            enddo
        enddo
        !$omp end do
    !$omp end parallel
    nEntries = outputIdx
end subroutine LKTomGrid