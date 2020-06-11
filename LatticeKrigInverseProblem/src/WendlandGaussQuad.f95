module QuadFunctions
    implicit none
    private
    public :: IntervalGaussQuad, WendlandEval

    contains
        !conducts numerical integration on a Wendland function along a line
        !offset is the distance between the line and the center
        !a and b are the start and end points of integration
        !quadPoints and quadWeights are the points to evaluate the function and the relative weights to estimate the
        !integral
        !nPoints is the number of points in the quadrature rule
        !nIntervals is the number of intervals to split the integral into (more intervals means less error)
        function IntervalGaussQuad(offset, a, b, quadPoints, quadWeights, nPoints, nIntervals) result(integralValue)
            integer, intent(in) :: nPoints, nIntervals
            integer :: intervalIdx
            double precision, intent(in) :: offset, a, b, quadPoints(nPoints), quadWeights(nPoints)
            double precision :: intervalLength, evalPoints(nPoints, nIntervals), evalDists(nPoints, nIntervals)
            double precision :: integralValue

            intervalLength = (b - a) / nIntervals
            do intervalIdx = 1,nIntervals
                evalPoints(:,intervalIdx) = a + intervalLength*(quadPoints(:) + intervalIdx-1)
            enddo
            evalDists = sqrt(evalPoints(:,:)**2 + offset**2)
            evalDists = WendlandEval(evalDists, 5, nIntervals)
            integralValue = 0
            do intervalIdx = 1,nIntervals
                integralValue = integralValue + sum(quadWeights * evalDists(:,intervalIdx))
            enddo
            integralValue = integralValue * intervalLength
        end function IntervalGaussQuad

        !evaluate the C4 Wendland function at every point in points, which is an nRows x nCols matrix
        function WendlandEval(points, nRows, nCols) result(values)
            integer, intent(in) :: nRows, nCols
            double precision, intent(in) :: points(nRows, nCols)
            double precision :: values(nRows, nCols)

            values = 1.d0/3 * (35*points**2 + 18*points + 3) * (1-points)**6
        end function WendlandEval
end module QuadFunctions

subroutine WendlandGaussQuad(offsets, completions, nEntries, nLevels, entries)
    use QuadFunctions
    implicit none
    integer, intent(in) :: nEntries, nLevels
    integer :: idx, intervalIdx, nQuad

    double precision, intent(in) :: offsets(nEntries), completions(2, nEntries)
    double precision :: intervalLength, xBound(nEntries), a, b
    double precision, intent(out) :: entries(nEntries)

    !5-point gaussian quadrature
    double precision, parameter :: quadPoints(5) = (/ 0.5 - sqrt(5+2*sqrt(10.d0/7))/6, 0.5 - sqrt(5-2*sqrt(10.d0/7))/6, 0.5d0, &
            0.5 + sqrt(5-2*sqrt(10.d0/7))/6, 0.5 + sqrt(5+2*sqrt(10.d0/7))/6 /)
    double precision, parameter :: quadWeights(5) = (/ (322-13*sqrt(70.d0))/1800, (322+13*sqrt(70.d0))/1800, 64.d0/225, &
            (322+13*sqrt(70.d0))/1800, (322-13*sqrt(70.d0))/1800 /)

    nQuad = size(quadPoints)
    entries(:) = 0
    xBound = sqrt(1-offsets**2)
    !$omp parallel private(a, b)
    !$omp do
    do idx = 1, nEntries
        if (completions(1,idx) == completions(2,idx)) then
            cycle
        endif
        !a and b are the beginning/end points of the interval
        a = xBound(idx) * (2*completions(1,idx) - 1)
        b = xBound(idx) * (2*completions(2,idx) - 1)
        !split the integral around 0, since the wendland function isn't analytic there
        if (a < 0 .and. b > 0) then
            entries(idx) = IntervalGaussQuad(offsets(idx), a, 0.d0, quadPoints, quadWeights, nQuad, nLevels) &
                          +IntervalGaussQuad(offsets(idx), 0.d0, b, quadPoints, quadWeights, nQuad, nLevels)
        else
            entries(idx) = IntervalGaussQuad(offsets(idx), a, b, quadPoints, quadWeights, nQuad, nLevels)
        endif
    enddo
    !$omp end do
    !$omp end parallel
end subroutine WendlandGaussQuad


