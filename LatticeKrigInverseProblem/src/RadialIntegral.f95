module AcosIntTime
    implicit none
    real :: t
end module AcosIntTime


module RadialIntegralFunctions
    implicit none
    private
    public :: CrossProd, SignedAngle, TriangleIntegral, WendlandEvalR

    contains


    !computes the cross product of two 2D vectors, returning the result as a scalar
    function CrossProd(v1, v2) result(s)
        double precision, intent(in) :: v1(2), v2(2)
        double precision :: s
        s = v1(1) * v2(2) - v1(2) * v2(1)
    end function CrossProd


    !computes the signed angle from L1 to a point, with the origin as the vertex of the angle
    !the returned angle is positive if the segment from L1 to L2 goes counterclockwise
        function SignedAngle(L1, L2) result(angle)
            double precision, intent(in) :: L1(2), L2(2)
            double precision :: angle
            angle = asin(crossProd(L1, L2) / sqrt(sum(L1**2) * sum(L2**2)))
            if (sum(L1 * L2) < 0) then
                angle = sign(1.d0, angle) * (3.14159265358979323846d0 - abs(angle))
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


    !integrates the C4 Wendland function centered at the origin with radius 1 over the triangle with vertices at
    !the origin, E1, and E2
    function TriangleIntegral(E1, E2) result(integralValue)
        double precision, intent(in) :: E1(2), E2(2)
        double precision :: L(2), perpDist, closerDist, longerDist, maxDist, longerPoint(2), longerAngle, totalSignedAngle
        double precision :: integralValue
        !5-point gaussian quadrature
        double precision, parameter :: quadPoints(5) = (/0.046910077030668d0, 0.230765344947158d0, 0.5d0, &
                0.769234655052841d0, 0.953089922969332d0/)
        double precision, parameter :: quadWeights(5) = (/0.118463442528095d0, 0.239314335249683d0, &
                0.284444444444444d0, 0.239314335249683d0, 0.118463442528095d0/)

        L = E2 - E1
        perpDist = abs(CrossProd(E1, L)) / sqrt(sum(L**2))
        closerDist = sqrt(min(sum(E1**2), sum(E2**2)))
        longerDist = sqrt(max(sum(E1**2), sum(E2**2)))
        maxDist = min(1d0, longerDist)
        if (closerDist < perpDist) then
            closerDist = perpDist
        endif
        if (longerDist <= closerDist) then
            integralValue = 0d0
            return
        endif
        if(sum(E1**2) > sum(E2**2)) then
            longerPoint = E1
        else
            longerPoint = E2
        endif
        longerAngle = acos(perpDist / longerDist)
        totalSignedAngle = SignedAngle(E1, E2)
        !integrate to the closer endpoint
        integralValue = totalSignedAngle * closerDist * sum(quadWeights * WendlandEvalR(quadPoints * closerDist, 5))
        !integrate for the radii between the two endpoints
        integralValue = integralValue + sign(perpDist*acosIntegral(closerDist/perpDist, maxDist/perpDist, perpDist,&
                longerAngle), totalSignedAngle)
    end function TriangleIntegral


    !split the integral into several different intervals: integrate from 1 to 2 using a transformation,
    !then integrate from 2 to 10, 10 to 100, and 100 to infinity directly. All integrals use 11-point
    !Gaussian quadrature
    function acosIntegral(xMin, xMax, perpDist, longerAngle) result(integralValue)
        double precision, intent(in) :: xMin, xMax, perpDist, longerAngle
        double precision :: intervalMin, intervalMax, tMin, tMax, tPoints(11), xPoints(11)
        double precision :: integralValue
        !11-point Gaussian quadrature
        double precision, parameter :: quadPoints(11) = (/ 0.010885670926972d0, 0.056468700115952d0, 0.134923997212975d0,&
                0.240451935396594d0, 0.365228422023827d0, 0.500000000000000d0, 0.634771577976172d0, 0.759548064603406d0,&
                0.865076002787025d0, 0.943531299884048d0, 0.989114329073028d0 /)
        double precision, parameter :: quadWeights(11) = (/ 0.027834283558087d0, 0.062790184732452d0, 0.093145105463867d0,&
                0.116596882295995d0, 0.131402272255123d0, 0.136462543388950d0, 0.131402272255123d0, 0.116596882295995d0,&
                0.093145105463867d0, 0.062790184732452d0, 0.027834283558087d0 /)

        integralValue = 0d0
        !integrate over the intersection of the intervals [1, 2] and [xMin, xMax]
        if (xMin <= 2) then
            intervalMax = min(2d0, xMax)
            tMin = acos(1/xMin)
            tMax = acos(1/intervalMax)
            tPoints = tMin + (tMax - tMin) * quadPoints
            xPoints = 1/cos(tPoints)
            integralValue = (tMax - tMin) * &
                    sum(quadWeights * WendlandEvalR(perpDist*xPoints, 11) * (longerAngle - tPoints) * sin(tPoints) * xPoints**2)
            if (xMax <= 2) then !the integral fits entirely in the interval [1, 2] so we're done
                return
            endif
        endif
        !integrate over the intersection of the intervals [2, 10] and [xMin, xMax]
        if (xMin <= 10) then
            intervalMin = max( 2d0, xMin)
            intervalMax = min(10d0, xMax)
            xPoints = intervalMin + (intervalMax - intervalMin) * quadPoints
            integralValue = integralValue + (intervalMax - intervalMin) * &
                    sum(quadWeights * WendlandEvalR(perpDist*xPoints, 11) * (longerAngle - acos(1/xPoints)))
            if (xMax <= 10) then
                return
            endif
        endif
        !integrate over the intersection of the intervals [10, 100] and [xMin, xMax]
        if (xMin <= 100) then
            intervalMin = max( 10d0, xMin)
            intervalMax = min(100d0, xMax)
            xPoints = intervalMin + (intervalMax - intervalMin) * quadPoints
            integralValue = integralValue + (intervalMax - intervalMin) * &
                    sum(quadWeights * WendlandEvalR(perpDist*xPoints, 11) * (longerAngle - acos(1/xPoints)))
            if (xMax <= 100) then
                return
            endif
        endif

        !integrate over the intersection of the intervals [100, +Inf] and [xMin, xMax]
        intervalMin = 100d0
        intervalMax = xMax
        xPoints = intervalMin + (intervalMax - intervalMin) * quadPoints
        integralValue = integralValue + (intervalMax - intervalMin) * &
                sum(quadWeights * WendlandEvalR(perpDist*xPoints, 11) * (longerAngle - acos(1/xPoints)))
    end function acosIntegral

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

        values = 1d0/3 * (35d0*points**2 + 18d0*points + 3d0) * points * (1d0-points)**6
        !replace wendland with 1 for debugging
        !values = points
    end function WendlandEvalR
end module RadialIntegralFunctions



subroutine RadialIntegral(vertices, nVerts, points, nPoints, ranges, rangeReps, nRanges, entries)
    use RadialIntegralFunctions
    implicit none
    integer, intent(in) :: nVerts, nPoints, nRanges, rangeReps(nRanges)
    integer :: rangeStarts(nRanges), rangeIdx, pointIdx, vertIdx

    double precision, intent(in) :: vertices(2, nVerts), points(2, nPoints), ranges(nRanges)
    double precision :: fullIntegral, range, range2, P(2), L1(2), L2(2), LP(2), L(2), t, closestDist2, integralValue
    double precision, intent(out) :: entries(nPoints)

    !for each point
    !for each edge
    !find the weight function w(r) for the triangle defined by the point and edge

    !fullIntegral is the integral from 0 to 1 of r times the wendland function
    fullIntegral = 1d0/18
    !when the wendland function is replaced by 1 for debugging
    !fullIntegral = 1d0/2
    rangeStarts(1) = 0
    do rangeIdx = 1, (nRanges-1)
        rangeStarts(rangeIdx+1) = rangeStarts(rangeIdx) + rangeReps(rangeIdx)
    enddo
    do rangeIdx = 1, nRanges
        range = ranges(rangeIdx)
        range2 = range**2
        !$omp parallel private(P, L1, L2, LP, L, t, integralValue, closestDist2)
        !$omp do
        do pointIdx = (rangeStarts(rangeIdx) + 1), (rangeStarts(rangeIdx) + rangeReps(rangeIdx))
            P = points(:,pointIdx)
            integralValue = 0
            do vertIdx = 1, nVerts
                L1 = vertices(:,vertIdx)
                if (vertIdx == nVerts) then
                    L2 = vertices(:,1)
                else
                    L2 = vertices(:,vertIdx+1)
                endif

                !check to make sure the point and edge aren't collinear (which would lead to floating point errors later)
                !recall that for vectors u and v, ||u x v|| = ||u|| ||v|| sin(a), where a is the angle between u and v
                !this check makes the code run about 50% slower, so i'm taking it out
                !if (CrossProd(L1-P, L2-P)**2 < 0.0000000001 * sum((L1-P)**2) * sum((L2-P)**2)) then
                !    cycle
                !endif
                L = L2 - L1
                t = sum((P - L1) * L) / sum(L**2)
                t = max(0d0, t)
                t = min(1d0, t)
                closestDist2 = sum((L1 + L*t - P)**2)
                if (closestDist2 > range2) then
                    integralValue = integralValue + SignedAngle(L1-P, L2-P) * fullIntegral * range2
                else
                    if (0 < t .and. t < 1) then
                        !break into two triangles by adding a vertex on the part of the line segment closest to the point
                        LP = L1 + t*L
                        integralValue = integralValue + TriangleIntegral((L1-P)/range, (LP-P)/range) * range2 &
                                                      + TriangleIntegral((LP-P)/range, (L2-P)/range) * range2
                    else
                        integralValue = integralValue + TriangleIntegral((L1-P)/range, (L2-P)/range) * range2
                    endif
                endif
            enddo
            !this is a bit sloppy; the more correct way would be to flip the sign of integralvalue iff the polygon is
            !defined clockwise, which we could test by finding a point inside the polygon and adding up totalSignedAngle
            !for every edge around that point; if it's -2*pi the polygon is defined CW, if it's 2*pi the polygon is CCW
            !however, as long as the basis function is nonnegative (like the Wendland is) then this works just as well
            entries(pointIdx) = abs(integralValue)
        enddo
        !$omp end do
        !$omp end parallel
    enddo
end subroutine RadialIntegral