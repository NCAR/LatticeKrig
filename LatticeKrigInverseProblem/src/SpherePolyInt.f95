module SpherePolyFunctions
    implicit none
    private
    public :: Normalize, CrossProd, SphereInterp, PointsInPolygon, SegmentIntegral, AdaptiveSegmentIntegral, DivergenceIntegral, &
            AdaptiveDivergenceIntegral, IntWendlandEval, ComputeFullIntegral, WendlandEval
    contains

    !normalize a vector, i.e. divide by its 2-norm
    subroutine Normalize(v)
        double precision, intent(inout) :: v(3)
        v = v / sqrt(sum(v**2))
    end subroutine Normalize


    !3D vector cross product
    function CrossProd(v1, v2) result(v3)
        double precision, intent(in) :: v1(3), v2(3)
        double precision :: v3(3)
        v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
        v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
        v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
    end function CrossProd


    !Linear interpolation over the surface of a sphere
    function SphereInterp(L1, L1Perp, a, b, t, nPoints) result(interpPoints)
        integer, intent(in) :: nPoints
        integer :: tIdx
        double precision, intent(in) :: L1(3), L1Perp(3), a, b, t(nPoints)
        double precision :: segmentLength
        double precision :: interpPoints(3, nPoints)
        segmentLength = b-a
        do tIdx = 1, nPoints
            interpPoints(:,tIdx) = L1*cos(a+t(tIdx)*segmentLength) + L1Perp*sin(a+t(tIdx)*segmentLength)
        enddo
    end function SphereInterp


    !determines whether each point in points is inside of the polygon on the surface of a sphere with the given vertices
    !we do this by requiring a point X that is known to be outside the polygon, then counting the number of times the arc
    !from X to each point intersects the polygon. The point is inside the polygon if and only if this number is odd
    !the returned value is an array of integers, where the nth entry is 1 if the nth point is inside the polygon, 0 otherwise.
    function PointsInPolygon(vertices, nVerts, points, nPoints, X) result(inPolygon)
        integer, intent(in) :: nVerts, nPoints
        integer :: pointIdx, vertIdx
        integer :: inPolygon(nPoints)
        double precision, intent(in) :: vertices(3, nVerts), points(3, nPoints), X(3)
        double precision :: P(3), v1(3), v2(3), v2Proj(nVerts), L1(3), L2(3), y1, y2, ixn(3)
        inPolygon(:) = 0
        do pointIdx = 1, nPoints
            !creating orthonormal vectors v1, v2, and P such that X is in the span of P and v1, so v2 is orthogonal to the
            !plane containing X, P, and the origin. we will project vertices onto v2 to determine which side of the arc from
            !X to P they are on
            P = points(:, pointIdx)
            if (sum(P*X) >= 1) cycle !P and X are the same point, and X is outside the polygon, so P must be too
            v1 = X - sum(P*X) * P
            call Normalize(v1)
            v2 = CrossProd(P, v1)
            v2Proj = matmul(v2, vertices)
            do vertIdx = 1, nVerts
                L1 = vertices(:, vertIdx)
                y1 = v2Proj(vertIdx)
                if (vertIdx == nVerts) then
                    L2 = vertices(:, 1)
                    y2 = v2Proj(1)
                else
                    L2 = vertices(:, vertIdx+1)
                    y2 = v2Proj(vertIdx+1)
                endif
                !both endpoints of the polygon's edge are on the same side of the arc from X to P, so
                !the arcs must not intersect
                if (y1 * y2 > 0d0) cycle
                !edge case handling when a vertex of the polygon is exactly on the arc from X to P; treat these points
                !as being to the right of the arc
                if (y1 == 0d0) then
                    if (y2 >= 0d0) cycle !both points are to the right of the arc
                    ixn = y1
                elseif (y2 == 0d0) then
                    if (y1 >= 0d0) cycle !both points are to the right of the arc
                    ixn = y2
                else
                    !compute the intersection of the polygon edge with the great circle defined by X and P
                    ixn = (y2*L1 - y1*L2) * sign(1d0, y2)
                    call Normalize(ixn)
                endif
                !if the intersection point is closer to P than X is, and it's in the this intersection point is on the arc from X to P
                !so we update the counter
                if (sum(P * ixn) > sum(P * X) .and. sum(v1 * ixn) > 0d0) inPolygon(pointIdx) = 1 - inPolygon(pointIdx)
            enddo
        enddo
    end function PointsInPolygon


    function SegmentIntegral(L1, L2, P, range) result(integralValue)
        double precision, intent(in) :: L1(3), L2(3), P(3), range
        double precision :: L1Perp(3), segmentLength, t, t2, LP(3)
        double precision :: integralValue
        double precision, parameter :: pi = acos(-1d0)
        !split the integral at the point closest to the center
        L1Perp = L2 - sum(L1*L2) * L1
        call Normalize(L1Perp)
        segmentLength = acos(sum(L1 * L2))
        t = atan2(sum(L1Perp*P), sum(L1*P))
        t2 = t - pi*sign(1d0, t)
        t = max(0d0, t)
        t = min(segmentLength, t) !clamp t to that so it stays on the curve
        t2 = max(0d0, t2)
        t2 = min(segmentLength, t2)
        if (0 < t .and. t < segmentLength) then
            LP = L1*cos(t) + L1Perp*sin(t)
            integralValue = DivergenceIntegral(L1, LP, P, range) + DivergenceIntegral(LP, L2, P, range)
        elseif (0 < t2 .and. t2 < segmentLength) then
            LP = L1*cos(t2) + L1Perp*sin(t2)
            integralValue = DivergenceIntegral(L1, LP, P, range) + DivergenceIntegral(LP, L2, P, range)
        else
            integralValue = DivergenceIntegral(L1, L2, P, range)
        endif
    end function SegmentIntegral


    function AdaptiveSegmentIntegral(L1, L2, P, range) result(integralValue)
        double precision, intent(in) :: L1(3), L2(3), P(3), range
        double precision :: L1Perp(3), segmentLength, t, t2
        double precision :: integralValue
        double precision, parameter :: pi = acos(-1d0)
        !split the integral at the points closest and farthest from the center
        L1Perp = L2 - sum(L1*L2) * L1
        call Normalize(L1Perp)
        segmentLength = acos(sum(L1 * L2))
        t = atan2(sum(L1Perp*P), sum(L1*P))
        t2 = t - pi*sign(1d0, t)
        t = max(0d0, t)
        t = min(segmentLength, t) !clamp t and t2 so they stay on the curve
        t2 = max(0d0, t2)
        t2 = min(segmentLength, t2)
        if (0 < t .and. t < segmentLength) then
            integralValue = AdaptiveDivergenceIntegral(L1, L1Perp, P, range, 0d0, t)&
                          + AdaptiveDivergenceIntegral(L1, L1Perp, P, range, t, segmentLength)
        elseif (0 < t2 .and. t2 < segmentLength) then
            integralValue = AdaptiveDivergenceIntegral(L1, L1Perp, P, range, 0d0, t2)&
                          + AdaptiveDivergenceIntegral(L1, L1Perp, P, range, t2, segmentLength)
        else
            integralValue = AdaptiveDivergenceIntegral(L1, L1Perp, P, range, 0d0, segmentLength)
        endif
    end function AdaptiveSegmentIntegral


    function DivergenceIntegral(L1, L2, P, range) result(integralValue)
        integer :: pointIdx
        double precision, intent(in) :: L1(3), L2(3), P(3), range
        double precision :: segmentLength, L1Perp(3), segmentPerp(3), integralPoints(3, 15), integralDists(15),&
                radialVectors(3, 15), integrandValues(15)
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
        segmentLength = acos(sum(L1 * L2))
        segmentPerp = -CrossProd(L1, L2)
        call Normalize(segmentPerp)
        L1Perp = L2 - sum(L2 * L1) * L1
        call Normalize(L1Perp)
        do pointIdx = 1,15
            integralPoints(:,pointIdx) = L1 * cos(quadPoints(pointIdx)*segmentLength) + &
                                     L1Perp * sin(quadPoints(pointIdx)*segmentLength)
            call Normalize(integralPoints(:,pointIdx))
            integralDists(pointIdx) = acos(sum(P * integralPoints(:,pointIdx)))
            radialVectors(:,pointIdx) = -(P - sum(P * integralPoints(:,pointIdx)) * integralPoints(:,pointIdx))
            call Normalize(radialVectors(:,pointIdx))
        enddo
        integrandValues = matmul(segmentPerp, radialVectors) * IntWendlandEval(range, integralDists, 15)
        integralValue = segmentLength * sum(quadWeightsKron * integrandValues)
    end function DivergenceIntegral


    !vectorized adaptive quadrature of the function we use in the divergence theorem over the interval from a to b
    function AdaptiveDivergenceIntegral(L1, L1Perp, P, range, a, b) result(integralValue)
        integer, parameter :: MAX_INTERVALS = 16
        integer :: nIntervals, newNIntervals, nPoints, pointIdx, intervalIdx, listIdx
        double precision, intent(in) :: L1(3), L1Perp(3), P(3), range, a, b
        double precision :: segmentPerp(3), endpoints(2, MAX_INTERVALS), segmentLength, abscissae(15), &
                integralPoints(3, 15*MAX_INTERVALS), integralDists(15*MAX_INTERVALS), intervalValues(MAX_INTERVALS), &
                intervalErrors(MAX_INTERVALS), radialVectors(3, 15*MAX_INTERVALS), integrandValues(15*MAX_INTERVALS), &
                integrandValueMatrix(15, MAX_INTERVALS), newEndpoints(2, MAX_INTERVALS), midpoint
        double precision :: integralValue
        double precision, parameter :: absTol = 1e-8, relTol = 1e-8
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

        segmentPerp = -CrossProd(L1, L1Perp)
        endpoints(:,1) = (/ a, b /)
        segmentLength = 2*(b-a)
        nIntervals = 1
        !this while loop runs the adaptive quadrature until the integral over every subinterval has been approximated
        !to within the absolute or relative error tolerances set above
        do while (nIntervals > 0)
            nPoints = 15*nIntervals
            segmentLength = segmentLength / 2
            do intervalIdx = 1, nIntervals
                abscissae = quadPoints * segmentLength + endpoints(intervalIdx, 1)
                do pointIdx = 1,15
                    listIdx = 15*(intervalIdx-1)+pointIdx
                    integralPoints(:, listIdx) = L1 * cos(abscissae(pointIdx)) + L1Perp * sin(abscissae(pointIdx));
                    radialVectors(:, listIdx) = -(P - sum(P * integralPoints(:,listIdx)) * integralPoints(:,listIdx))
                    call Normalize(radialVectors(:, listIdx))
                enddo
            enddo
            !evaluate the integrand function at each point, scale by the dot product of the line's normal
            !vector with the direction of the vector field (from Divergence Theorem)
            integralDists(1:nPoints) = acos(matmul(P, integralPoints(:, 1:nPoints)))
            integrandValues(1:nPoints) = matmul(segmentPerp, radialVectors(:, 1:nPoints)) * &
                    IntWendlandEval(range, integralDists(1:nPoints), nPoints)
            integrandValueMatrix = reshape(integrandValues, (/ 15, MAX_INTERVALS /))
            intervalValues = segmentLength * matmul(quadWeightsKron,  integrandValueMatrix)
            intervalErrors = abs(intervalValues - segmentLength * &
                    matmul(quadWeightsGauss, integrandValueMatrix))
            newNIntervals = 0
            do intervalIdx = 1, nIntervals
                !accept the approximation over the current interval if the interval is too small, the next interval buffer
                !is full, or the approximation meets an error tolerance
                if (segmentLength * 1e13 < (b-a) .or. newNIntervals >= MAX_INTERVALS .or. &
                        intervalErrors(intervalIdx) < absTol * segmentLength/(b-a) .or. &
                        intervalErrors(intervalIdx) < relTol * abs(intervalValues(intervalIdx))) then

                    integralValue = integralValue + intervalValues(intervalIdx)
                else
                    !split the current interval in half and evaluate those two intervals in the next step
                    midpoint = sum(endpoints(:, intervalIdx))/2
                    newEndpoints(1, newNIntervals+1) = endpoints(1, intervalIdx)
                    newEndpoints(2, newNIntervals+1) = midpoint
                    newEndpoints(1, newNIntervals+2) = midpoint
                    newEndpoints(2, newNIntervals+2) = endpoints(2, intervalIdx)
                    newNIntervals = newNIntervals + 2
                endif
            enddo
            nIntervals = newNIntervals
            endpoints = newEndpoints
        enddo
    end function AdaptiveDivergenceIntegral


    !evaluates the function we need for the Divergence Theorem: 1/sin(theta) * integral from 0 to theta of sin(t) *
    !W(t/radius) dt, where W is the C4 Wendland function below, and theta is an element of the array points
    function IntWendlandEval(radius, points, nPoints) result(values)
        integer, intent(in) :: nPoints
        integer :: pointIdx
        double precision, intent(in) :: radius, points(nPoints)
        double precision :: integralPoints(11)
        double precision :: values(nPoints)
        !11-point Gaussian quadrature, since the integrand is analytic
        double precision, parameter :: quadWeights(11) = (/ 0.027834283558087d0, 0.062790184732452d0, 0.093145105463867d0,&
                0.116596882295995d0, 0.131402272255123d0, 0.136462543388950d0, 0.131402272255123d0, 0.116596882295995d0,&
                0.093145105463867d0, 0.062790184732452d0, 0.027834283558087d0 /)
        double precision, parameter :: quadPoints(11) = (/ 0.010885670926972d0, 0.056468700115952d0, 0.134923997212975d0,&
                0.240451935396594d0, 0.365228422023827d0, 0.500000000000000d0, 0.634771577976172d0, 0.759548064603406d0,&
                0.865076002787025d0, 0.943531299884048d0, 0.989114329073028d0 /)
        do pointIdx = 1,nPoints
            !don't integrate outside the Wendland function
            integralPoints = quadPoints * min(points(pointIdx), radius)
            values(pointIdx) = min(points(pointIdx), radius) * &
                    sum(quadWeights * sin(integralPoints) * WendlandEval(integralPoints / radius, 11))
        enddo
        values(:) = values(:) / sin(points)
    end function IntWendlandEval


    !integrate the basis function with the given range over the whole sphere
    function ComputeFullIntegral(range) result(fullSphereIntegral)
        double precision, intent(in) :: range
        double precision :: abscissae(11)
        double precision :: fullSphereIntegral
        double precision, parameter :: pi = acos(-1d0)
        !11-point Gaussian quadrature, since the integrand is analytic
        double precision, parameter :: quadWeights(11) = (/ 0.027834283558087d0, 0.062790184732452d0, 0.093145105463867d0,&
                0.116596882295995d0, 0.131402272255123d0, 0.136462543388950d0, 0.131402272255123d0, 0.116596882295995d0,&
                0.093145105463867d0, 0.062790184732452d0, 0.027834283558087d0 /)
        double precision, parameter :: quadPoints(11) = (/ 0.010885670926972d0, 0.056468700115952d0, 0.134923997212975d0,&
                0.240451935396594d0, 0.365228422023827d0, 0.500000000000000d0, 0.634771577976172d0, 0.759548064603406d0,&
                0.865076002787025d0, 0.943531299884048d0, 0.989114329073028d0 /)
        abscissae = quadPoints * min(2*pi, range)
        fullSphereIntegral = 2*pi * min(2*pi, range) * sum(quadWeights * sin(abscissae) * WendlandEval(abscissae / range, 11))
    end function ComputeFullIntegral


    !evaluates the C4 Wendland function at each element of points; this function doesn't check that the entries in points
    !are less than 1, and will give incorrect results for entries greater than 1, because it's more efficient to guarantee
    !that elsewhere
    function WendlandEval(points, nPoints) result(values)
        integer, intent(in) :: nPoints
        double precision, intent(in) :: points(nPoints)
        double precision :: values(nPoints)
        values = 1.d0/3 * (35*points**2 + 18*points + 3) * (1-points)**6
        !replace Wendland with 1 for debugging
        !values(:) = 1
    end function WendlandEval
end module SpherePolyFunctions



subroutine SpherePolyInt(vertices, nVerts, points, nPoints, ranges, rangeReps, nRanges, entries)
    use SpherePolyFunctions
    implicit none
    integer, intent(in) :: nVerts, nPoints, nRanges, rangeReps(nRanges)
    integer :: rangeStarts(nRanges), rangeIdx, pointIdx, vertIdx, oppositeInPolyMask(nPoints)
    double precision, intent(in) :: vertices(3, nVerts), points(3, nPoints), ranges(nRanges)
    double precision :: extPoint(3), L1(3), L2(3), fullSphereIntegral
    double precision, intent(out) :: entries(nPoints)
    extPoint = -sum(vertices, 2)
    call Normalize(extPoint)
    oppositeInPolyMask = pointsInPolygon(vertices, nVerts, -1 * points, nPoints, extPoint)
    rangeStarts(1) = 0
    do rangeIdx = 1, (nRanges-1)
        rangeStarts(rangeIdx+1) = rangeStarts(rangeIdx) + rangeReps(rangeIdx)
    enddo
    do rangeIdx = 1, nRanges
        fullSphereIntegral = ComputeFullIntegral(ranges(rangeIdx))
        !$omp parallel private()
        !$omp do
        do pointIdx = (rangeStarts(rangeIdx) + 1), (rangeStarts(rangeIdx) + rangeReps(rangeIdx))
            do vertIdx = 1,nVerts
                L1 = vertices(:,vertIdx)
                if (vertIdx == nVerts) then
                    L2 = vertices(:,1)
                else
                    L2 = vertices(:,vertIdx+1)
                endif
                !entries(pointIdx) = entries(pointIdx) + SegmentIntegral(L1, L2, points(:,pointIdx), ranges(rangeIdx))
                entries(pointIdx) = entries(pointIdx) + AdaptiveSegmentIntegral(L1, L2, points(:,pointIdx), ranges(rangeIdx))
            enddo
            !if the point opposite the center is inside the polygon, the result computed above will be the negative integral
            !of the function over the exterior of the polygon, so we add the integral over the whole sphere to get the right
            !answer
            entries(pointIdx) = entries(pointIdx) + fullSphereIntegral * oppositeInPolyMask(pointIdx)
        enddo
        !$omp end do
        !$omp end parallel
    enddo
end subroutine SpherePolyInt
