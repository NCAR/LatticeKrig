module RevisedIntegralFunctions
    implicit none
    private
    public :: CrossProd, SignedAngle, SignedAngleVector, TriangleWeights, TriangleWeightsEvalGrid, WendlandEvalR

contains

    !quick approximation of acos(x), since that function was taking 75% of the execution time
    function qacos(x) result (y)
        double precision, intent(in) :: x
        double precision :: y
        !y = sqrt(2-2*x) + 0.1566*(1-x)**2 - 4/125 * x * (x-1)
        y = sqrt(2-2*x)
        y = y + (1-0.5*y**2 - x) / (y - 1d0/6*y**3)
        y = y + (1-0.5*y**2 - x) / (y - 1d0/6*y**3)
        !y = y + (cos(y) - x) / sin(y)
    end function qacos
    !computes the cross product of two 2D vectors, returning the result as a scalar
    function CrossProd(v1, v2) result(s)
        double precision, intent(in) :: v1(2), v2(2)
        double precision :: s

        s = v1(1) * v2(2) - v1(2) * v2(1)
    end function CrossProd


    function SignedAngle(L1, L2) result(angle)
        double precision, intent(in) :: L1(2), L2(2)
        double precision :: angle
        angle = asin(crossProd(L1, L2) / sqrt(sum(L1**2) * sum(L2**2)))
        if (sum(L1 * L2) < 0) then
            !atan2(0, -1) = pi
            angle = sign(1.d0, angle) * (atan2(0.d0, -1.d0) - abs(angle))
        endif
    end function SignedAngle


    !computes the signed angle from L1 to each element of points, with the origin as the vertex of the angle
    !the returned angle is positive if it goes counterclockwise
    function SignedAngleVector(L1, points, nPoints) result(angles)
        integer, intent(in) :: nPoints
        integer :: angleIdx
        double precision, intent(in) :: L1(2), points(2,nPoints)
        double precision :: angles(nPoints)
        do angleIdx = 1,nPoints
            angles(angleIdx) = SignedAngle(L1, points(:,angleIdx))
        enddo
    end function SignedAngleVector


    !evaluates the weight function over the triangle with vertices at the origin, E1, and E2
    function TriangleWeights(E1, E2, evalRanges, nRanges) result(weights)
        integer, intent(in) :: nRanges
        integer :: rangeIdx
        double precision, intent(in) :: E1(2), E2(2), evalRanges(nRanges)
        double precision :: L(2), perpDist, closerDist, longerDist, longerPoint(2), totalSignedAngle, maxAngle
        double precision :: weights(nRanges)

        L = E2 - E1
        totalSignedAngle = SignedAngle(E1, E2)
        closerDist = sqrt(min(sum(E1**2), sum(E2**2)))
        longerDist = sqrt(max(sum(E1**2), sum(E2**2)))
        !we assume later that perpDist <= closerDist, which would always be true except for floating point arithmetic, so
        !this min fixes those potential problems
        perpDist = min(closerDist, abs(CrossProd(E1, L)) / sqrt(sum(L**2)))
        maxAngle = acos(perpDist / longerDist)
        if(sum(E1**2) > sum(E2**2)) then
            longerPoint = E1
        else
            longerPoint = E2
        endif
        weights(:) = 0
        do rangeIdx = 1, nRanges
            if (evalRanges(rangeIdx) <= closerDist) then
                weights(rangeIdx) = totalSignedAngle
            elseif (evalRanges(rangeIdx) >= longerDist) then
                exit
            else
                !weights(rangeIdx) = sign(maxAngle - acos(perpDist / evalRanges(rangeIdx)), totalSignedAngle)
                weights(rangeIdx) = sign(maxAngle - qacos(perpDist / evalRanges(rangeIdx)), totalSignedAngle)
            endif
        enddo
    end function TriangleWeights


    function TriangleWeightsEvalGrid(E1, E2, evalRangeGap, nRanges) result(weights)
        integer, intent(in) :: nRanges
        integer :: startIdx, stopIdx, midIdx
        double precision, intent(in) :: E1(2), E2(2), evalRangeGap
        double precision :: L(2), perpDist, closerDist, longerDist, totalSignedAngle, maxAngle
        double precision :: weights(nRanges)

        L = E2 - E1
        totalSignedAngle = SignedAngle(E1, E2)
        closerDist = sqrt(min(sum(E1**2), sum(E2**2)))
        longerDist = sqrt(max(sum(E1**2), sum(E2**2)))
        !we assume later that perpDist <= closerDist, which would always be true except for floating point arithmetic, so
        !this min fixes those potential problems
        perpDist = min(closerDist, abs(CrossProd(E1, L)) / sqrt(sum(L**2)))
        maxAngle = acos(perpDist / longerDist)
        startIdx = floor(closerDist / evalRangeGap)
        stopIdx  = floor(longerDist / evalRangeGap)
        startIdx = max(startIdx, 0)
        stopIdx  = min(stopIdx, nRanges)
        if (startIdx > 0) weights(:startIdx) = 1
        if (stopIdx < nRanges) weights((stopIdx+1):) = 0
        do midIdx = (startIdx+1), stopIdx
            weights(midIdx) = sign(maxAngle - acos(perpDist / ((midIdx-1) * evalRangeGap)), totalSignedAngle)
        enddo
    end function TriangleWeightsEvalGrid



    function PointOnSegmentAtRadius(E1, E2, radius) result (P)
        double precision, intent(in) :: E1(2), E2(2), radius
        double precision :: L(2), tMid, t, P(2)

        !this function involves solving a quadratic equation; we are only interested in the solution between 0 and 1, and
        !the function will only be called in situatitions where exactly 1 solution in this interval is guaranteed to exist
        L = E2 - E1
        tMid = -sum(E1*L) / sum(L**2)
        t = tMid - sign(sqrt(radius**2 * sum(L**2) - (L(1) * E1(2) - L(2) * E1(1))**2), tMid - 0.5d0) / sum(L**2)
        P = E1 + L*t
    end function PointOnSegmentAtRadius


    !evaluate r times the C4 Wendland function at every point in points, which is a vector of length nEntries
    function WendlandEvalR(points, nEntries) result(values)
        integer, intent(in) :: nEntries
        double precision, intent(in) :: points(nEntries)
        double precision :: values(nEntries)

        values = 1.d0/3 * (35.d0*points**2 + 18.d0*points + 3.d0) * points * (1.d0-points)**6
        !replace wendland with 1 for debugging
        !values = points
    end function WendlandEvalR
end module RevisedIntegralFunctions



subroutine DeprecatedRadialIntegral(vertices, nVerts, points, nPoints, ranges, rangeReps, nRanges, entries, nEvals)
    use RevisedIntegralFunctions
    implicit none
    integer, intent(in) :: nVerts, nPoints, nRanges, rangeReps(nRanges), nEvals
    integer :: rangeStarts(nRanges), rangeIdx, pointIdx, vertIdx
    double precision, intent(in) :: vertices(2, nVerts), points(2, nPoints), ranges(nRanges)
    double precision :: quadWeights(nEvals), range, range2, P(2), L1(2), L2(2), LP(2), L(2), t, closestDist2
    double precision :: weightEvals(nEvals), evalRanges(nEvals), evalRangeGap
    double precision, intent(out) :: entries(nPoints)
    !double precision, intent(out) :: entries(nVerts, nPoints)

    evalRangeGap = 1.d0/(nEvals-1)
    evalRanges = (/ (rangeIdx, rangeIdx=0,(nEvals-1)) /) * evalRangeGap

    quadWeights(:) = 2
    quadWeights(1) = 1
    quadWeights(nEvals) = 1
    do rangeIdx = 2, nEvals-1, 2
        quadWeights(rangeIdx) = 4
    enddo

    entries(:) = -5.4321d0
    rangeStarts(1) = 0
    do rangeIdx = 1, (nRanges-1)
        rangeStarts(rangeIdx+1) = rangeStarts(rangeIdx) + rangeReps(rangeIdx)
    enddo
    do rangeIdx = 1, nRanges
        range = ranges(rangeIdx)
        range2 = range**2
        !$omp parallel private(P, L1, L2, LP, L, t, closestDist2, weightEvals)
        !$omp do
        do pointIdx = (rangeStarts(rangeIdx) + 1), (rangeStarts(rangeIdx) + rangeReps(rangeIdx))
            weightEvals(:) = 0
            P = points(:,pointIdx)
            do vertIdx = 1, nVerts
                L1 = vertices(:,vertIdx)
                if (vertIdx == nVerts) then
                    L2 = vertices(:,1)
                else
                    L2 = vertices(:,vertIdx+1)
                endif
                !finding the shortest distance from P to the line segment between L1 and L2
                L = L2 - L1
                t = dot_product(P - L1, L) / dot_product(L, L)
                t = max(0.d0, t)
                t = min(1.d0, t)
                closestDist2 = sum((L1 + t*L - P)**2)
                if (closestDist2 > range2) then
                    weightEvals = weightEvals + SignedAngle(L1-P, L2-P)
                    !weightEvals = SignedAngle(L1-P, L2-P)
                    !entries(vertIdx, pointIdx) = weightEvals(2)
                else
                    if (0.d0 < t .and. t < 1.d0) then
                        !break into two triangles by adding a vertex on the part of the line segment closest to the point
                        LP = L1 + t*L
                        weightEvals = weightEvals + TriangleWeights((L1-P)/range, (LP-P)/range, evalRanges, nEvals)&
                                + TriangleWeights((LP-P)/range, (L2-P)/range, evalRanges, nEvals)
                    else
                        weightEvals = weightEvals + TriangleWeights((L1-P)/range, (L2-P)/range, evalRanges, nEvals)
                    endif
                endif
            enddo
            entries(pointIdx) = sum(quadWeights * abs(weightEvals) * WendlandEvalR(evalRanges, nEvals)) * range2 * evalRangeGap/3
        enddo
        !$omp end do
        !$omp end parallel
    enddo
end subroutine DeprecatedRadialIntegral
