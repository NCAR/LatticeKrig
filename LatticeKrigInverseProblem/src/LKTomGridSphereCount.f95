subroutine LKTomGridSphereCount(points, nPoints, lines, nLines, ranges, rangeReps, nRanges, nEntries)
    implicit none
    integer, parameter :: dim = 3
    integer, intent(in) :: nPoints, nLines, nRanges, rangeReps(nRanges)
    integer :: rangeStarts(nRanges), outputIdx, rangeIdx, lineIdx, dimIdx, pointIdx
    integer, intent(out) :: nEntries

    double precision, intent(in) :: points(dim, nPoints), lines(2*dim, nLines), ranges(nRanges)
    double precision :: pointVec(dim), range, v1(dim), v2(dim), v1Perp(dim), v2Perp(dim), normal(dim)
    double precision :: dotProd, projRange

    !for each line
     !compute cross product of the line's endpoints to get unit normal vector
     !for each range
      !take the sine of the range to get the euclidean distance in the direction of the normal vector that guarantees
      !that a line misses a point
      !for each point
       !take the inner product of that point with the normal vector; if its absolute value is greater than the
       !sine of the current range, continue (the line must miss the point)
       !otherwise, check to see if the point misses the line segment:
       !if the arccos of the inner product of the point with either of the two line points is less than the range, the
       !point is in range of that endpoint; otherwise, the circle must either touch only the interior of the line or not
       !touch the line at all. We test this by projecting the point onto the vectors perpendicular to L1 and L2,
       !pointing inward in the plane of L1 and L2. if either projection is negative, the point must be closer to
       !an endpoint than any other point on the line segment, so it doesn't touch; otherwise, the point must be within
       !range of the line, and on the section of the sphere where its closer point to the line segment is on the
       !interior, so it must touch the line.

       !pseudocode:
        !check if circle is too far away from whole line
        !check if circle is close enough to an endpoint
        !check if circle is along the interior of the line segment
       !count up total number of points, report that as nEntries

    outputIdx = 0

    ! determining which points will have which ranges
    rangeStarts(1) = 0
    do rangeIdx = 1, (nRanges-1)
        rangeStarts(rangeIdx+1) = rangeStarts(rangeIdx) + rangeReps(rangeIdx)
    enddo

    !$omp parallel private(pointVec, v1, v2, dotProd, v1Perp, v2Perp, normal, range, projRange)
    !$omp do reduction(+:outputIdx)
    do lineIdx = 1, nLines
        v1 = lines(:dim, lineIdx)
        v2 = lines(dim+1:, lineIdx)
        dotProd = sum(v1 * v2)
        v1Perp = v2 - dotProd * v1
        v2Perp = v1 - dotProd * v2
        normal(1) = v1(2) * v2(3) - v1(3) * v2(2)
        normal(2) = v1(3) * v2(1) - v1(1) * v2(3)
        normal(3) = v1(1) * v2(2) - v1(2) * v2(1)
        normal = normal / sqrt(sum(normal * normal))
        do rangeIdx = 1, nRanges
            range = ranges(rangeIdx)
            projRange = sin(range)
            do pointIdx = (rangeStarts(rangeIdx)+1), (rangeStarts(rangeIdx) + rangeReps(rangeIdx))
                pointVec = points(:, pointIdx)
                !check 1
                ! the cosine check is to make sure that range is less than pi/2; if range is pi/2 or greater, then
                ! every possible line will intersect the circle at some point, so this if statement should always be false
                if (cos(range) > 0 .and. abs(sum(pointVec * normal)) > projRange) then
                    cycle
                endif

                !check 2
                if (acos(sum(v1 * pointVec)) < range .or. acos(sum(v2 * pointVec)) < range) then
                    outputIdx = outputIdx + 1
                    cycle
                endif

                !check 3
                if (sum(v1Perp * pointVec) > 0 .and. sum(v2Perp * pointVec) > 0) then
                    outputIdx = outputIdx + 1
                endif
            enddo
        enddo
    enddo
    !$omp end do
    !$omp end parallel

    nEntries = outputIdx
end subroutine LKTomGridSphereCount
