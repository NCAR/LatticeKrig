! identical to LKTomPoints3D except without the 3rd dimension on points, so there are only 2 do loops and we don't set
! dim3Start or points(3, outputIdx)

subroutine LKTomPoints2D(coordinates, levelSizes, nLevels, points, nPoints)
    implicit none
    integer, parameter :: dimension = 2
    integer, intent(in) :: nLevels, nPoints, levelSizes(dimension, nLevels)
    integer :: outputIdx, levelIdx, dimIdx, dim1Idx, dim2Idx, levelStart, dim1Start, dim2Start, &
            starts(dimension, nLevels)

    double precision, intent(in) :: coordinates(sum(levelSizes))
    double precision, intent(out) :: points(dimension, nPoints)

    outputIdx = 1
    do levelIdx = 1, nLevels
        levelStart = sum(levelSizes(:,1:(levelIdx-1)))
        do dimIdx = 1, dimension
            starts(dimIdx, levelIdx) = levelStart + sum(levelSizes(1:(dimIdx-1), levelIdx))
        enddo
    enddo
    do levelIdx = 1, nLevels
        dim1Start = starts(1, levelIdx)
        dim2Start = starts(2, levelIdx)
        do dim1Idx = 1, levelSizes(1, levelIdx)
            do dim2Idx = 1, levelSizes(2, levelIdx)
                points(1, outputIdx) = coordinates(dim1Start + dim1Idx)
                points(2, outputIdx) = coordinates(dim2Start + dim2Idx)
                outputIdx = outputIdx + 1
            enddo
        enddo
    enddo
end subroutine LKTomPoints2D