!takes the format of the lattice stored in an LKinfo object by LKrigSetup and creates the full list of all the lattice
!points that we need to measure distances against in LKTomGrid.

subroutine LKTomPoints3D(coordinates, levelSizes, nLevels, points, nPoints)
    implicit none
    integer, parameter :: dimension = 3
    integer, intent(in) :: nLevels, nPoints, levelSizes(dimension, nLevels)
    integer :: outputIdx, levelIdx, dimIdx, dim1Idx, dim2Idx, dim3Idx, levelStart, dim1Start, dim2Start, dim3Start, &
            starts(dimension, nLevels)

    double precision, intent(in) :: coordinates(sum(levelSizes))
    double precision, intent(out) :: points(dimension, nPoints)

    outputIdx = 1
    !computing which index of coordinates each dimension for each level starts at
    do levelIdx = 1, nLevels
        levelStart = sum(levelSizes(:,1:(levelIdx-1)))
        do dimIdx = 1, dimension
            starts(dimIdx, levelIdx) = levelStart + sum(levelSizes(1:(dimIdx-1), levelIdx))
        enddo
    enddo
        do levelIdx = 1, nLevels
            dim1Start = starts(1, levelIdx)
            dim2Start = starts(2, levelIdx)
            dim3Start = starts(3, levelIdx)
            do dim1Idx = 1, levelSizes(1, levelIdx)
                do dim2Idx = 1, levelSizes(2, levelIdx)
                    do dim3Idx = 1, levelSizes(3, levelIdx)
                        !fetch the coordinates and write them into the next point
                        points(1, outputIdx) = coordinates(dim1Start + dim1Idx)
                        points(2, outputIdx) = coordinates(dim2Start + dim2Idx)
                        points(3, outputIdx) = coordinates(dim3Start + dim3Idx)
                        outputIdx = outputIdx + 1
                    enddo
                enddo
            enddo
        enddo
end subroutine LKTomPoints3D