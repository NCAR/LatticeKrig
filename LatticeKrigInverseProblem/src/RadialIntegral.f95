module RadialIntegralFunctions
    implicit none
    private
    public :: CrossProd, SignedAngle, SegmentIntegral, DivergenceIntegral, IntRWendlandEval

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
            angle = sign(1.d0, angle) * (3.14159265358979323d0 - abs(angle))
        endif
    end function SignedAngle


    !sets up to integrate the C4 Wendland function centered at the origin with radius 1 over the triangle with vertices at
    !the origin, L1, and L2 by using the divergence theorem
    function SegmentIntegral(L1, L2) result(integralValue)
        double precision, intent(in) :: L1(2), L2(2)
        double precision :: L(2), splitPoint(2), closerDist, longerDist, angleSign
        double precision :: integralValue
        double precision, parameter :: fullIntegral = 1d0/18

        L = L2 - L1
        angleSign = CrossProd(L1, L2)
        closerDist = sqrt(min(sum(L1**2), sum(L2**2)))
        longerDist = sqrt(max(sum(L1**2), sum(L2**2)))
        !this should only happen if the two points L1 and L2 are on top of each other
        if (longerDist <= closerDist) then
            integralValue = 0d0
            return
        endif
        !if the line segment goes out of the basis function's support, split the segment into the parts inside and outside
        !the support
        !this assumes the shorterDist is <= 1, which needs to be guaranteed by the function calling this
        if (longerDist > 1d0) then
            splitPoint = PointOnSegmentAtRadius(L1, L2, 1d0)
            !use divergence theorem to integrate inside the basis function; use geometry to integrate outside
            !the functions DivergenceIntegral and signedAngle are antisymmetric, so we have to make sure we use
            !L1, splitPoint, and L2 in that order
            if(sum(L1**2) > sum(L2**2)) then
                integralValue = DivergenceIntegral(splitPoint, L2) + fullIntegral * signedAngle(L1, splitPoint)
            else
                integralValue = DivergenceIntegral(L1, splitPoint) + fullIntegral * signedAngle(splitPoint, L2)
            endif
        else
            integralValue = DivergenceIntegral(L1, L2)
        endif
    end function SegmentIntegral

    !computes the integral over the line segment from L1 to L2 of the function we get from the divergence theorem,
    !f(r) = 1/r * integral from 0 to r of t phi(t) dt, dotted with the unit normal vector pointing to the right of the
    !segment (outward if the polygon is specified CCW)
    function DivergenceIntegral(L1, L2) result(integralValue)
        integer :: pointIdx
        double precision, intent(in) :: L1(2), L2(2)
        double precision :: L(2), LPerp(2), integralPoints(2,15), integralDists(15)
        double precision :: integralValue
        !G7, K15 Gauss-Kronrod quadrature
        double precision, parameter :: quadPoints(15) = (/ 0.004272314439594d0, 0.025446043828621d0, 0.067567788320115d0,&
                0.129234407200303d0, 0.206956382266154d0, 0.297077424311301d0, 0.396107522496051d0, 0.500000000000000d0,&
                0.603892477503949d0, 0.702922575688699d0, 0.793043617733846d0, 0.870765592799697d0, 0.932432211679885d0,&
                0.974553956171379d0, 0.995727685560406d0 /)
        double precision, parameter :: quadWeightsKron(15) = (/ 0.011467661005265d0, 0.031546046314989d0, 0.052395005161125d0,&
                0.070326629857763d0, 0.084502363319634d0, 0.095175289032393d0, 0.102216470037649d0, 0.104741070542364d0,&
                0.102216470037649d0, 0.095175289032393d0, 0.084502363319634d0, 0.070326629857763d0, 0.052395005161125d0,&
                0.031546046314989d0, 0.011467661005265d0 /)
        double precision, parameter :: quadWeightsGauss(15) = (/ 0d0, 0.064742483084435d0, 0d0, 0.139852695744638d0, 0d0,&
                0.190915025252559d0, 0d0, 0.208979591836735d0, 0d0, 0.190915025252559d0, 0d0, 0.139852695744638d0, 0d0,&
                0.064742483084435d0, 0d0 /)

        L = L2 - L1
        LPerp(1) = L(2)
        LPerp(2) = -L(1)
        do pointIdx = 1, 15
            integralPoints(:,pointIdx) = L1 + L*quadPoints(pointIdx)
            integralDists(pointIdx) = sqrt(sum(integralPoints(:,pointIdx)**2))
        enddo
        !to save time (and since the integrand is analytic so the computed result is very accurate), we don't do
        !adaptive quadrature or error estimation
        integralValue = sum(quadWeightsKron * &
                IntRWendlandEval(integralDists, 15) * matmul(LPerp, integralPoints) / integralDists)
    end function DivergenceIntegral


    !finds a point on the line through L1 and L2 that's on the circle with the given radius centered at the origin
    !if no such point exists, a floating point error will occur; if two such points exist, the one closer to the midpoint
    !of L1 and L2 is returned
    function PointOnSegmentAtRadius(L1, L2, radius) result (P)
        double precision, intent(in) :: L1(2), L2(2), radius
        double precision :: L(2), tMid, t, P(2)

        !this function involves solving a quadratic equation; we are only interested in the solution between 0 and 1, and
        !the function will only be called in situatitions where exactly 1 solution in this interval is guaranteed to exist
        L = L2 - L1
        tMid = -sum(L1*L) / sum(L**2)
        t = tMid - sign(sqrt(radius**2 * sum(L**2) - (L(1) * L1(2) - L(2) * L1(1))**2), tMid - 0.5d0) / sum(L**2)
        P = L1 + L*t
    end function PointOnSegmentAtRadius


    !evaluates 1/R times the integral from 0 to R of t phi(t) dt for each R in x
    function IntRWendlandEval(x, nEntries) result(values)
        integer, intent(in) :: nEntries
        double precision, intent(in) :: x(nEntries)
        double precision :: values(nEntries)

        values = 1d0/18 * x * (9 - 42*x**2 + 210*x**4 - 384*x**5 + 315*x**6 - 128*x**7 + 21*x**8)
        !check to make sure we didn't call the function for a distance greater than 1
        !it's more efficient to prevent that elsewhere in the code than to do it here, so this is commented out
        !do entryIdx = 1, nEntries
        !    if(points(entryIdx) > 1d0) values(entryIdx) = 1d0/18 / points(entryIdx)
        !enddo

        !replace wendland with constant function for debugging
        !values = 0.5d0 * x
    end function IntRWendlandEval
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
                    integralValue = integralValue + SignedAngle(L1-P, L2-P) * fullIntegral
                else
                    LP = L1 + t*L
                    integralValue = integralValue + SegmentIntegral((L1-P)/range, (LP-P)/range)&
                                                  + SegmentIntegral((LP-P)/range, (L2-P)/range)
                endif
            enddo
            !this is a bit sloppy; the more correct way would be to flip the sign of integralvalue iff the polygon is
            !defined clockwise, which we could test by adding the signed exterior angles for every vertex:
            !if it's -2*pi the polygon is defined CW, if it's 2*pi the polygon is CCW
            !however, as long as the basis function is nonnegative (like the Wendland is) then this works just as well
            entries(pointIdx) = abs(integralValue) * range2
        enddo
        !$omp end do
        !$omp end parallel
    enddo
end subroutine RadialIntegral
