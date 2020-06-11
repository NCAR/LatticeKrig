module PointNearPolygonFunctions
    implicit none
    private
    public :: CrossProd, SignedAngle, MinSquaredDistance
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
            !acos(-1) = pi
            angle = sign(1.d0, angle) * (acos(-1d0) - abs(angle))
        endif
    end function SignedAngle


    !computes the minimum squared distance from the origin to the line segment with endpoints L1 and L2
    function MinSquaredDistance(L1, L2) result (dist2)
        double precision, intent(in) :: L1(2), L2(2)
        double precision :: L(2), t, dist2

        L = L2 - L1
        t = -sum(L1 * L) / sum(L * L)
        t = max(0d0, t)
        t = min(1d0, t)
        dist2 = sum((L1 + L*t)**2)
    end function MinSquaredDistance
end module PointNearPolygonFunctions

subroutine PointNearPolygon(vertices, nVerts, points, nPoints, ranges, rangeReps, nRanges, entries)
    use PointNearPolygonFunctions
    implicit none
    integer, intent(in) :: nVerts, nPoints, nRanges, rangeReps(nRanges)
    integer :: rangeStarts(nRanges), rangeIdx, pointIdx, vertIdx
    integer, intent(out) :: entries(nPoints)

    double precision, intent(in) :: vertices(2, nVerts), points(2, nPoints), ranges(nRanges)
    double precision :: range, range2, P(2), L1(2), L2(2), sweptAngle
    double precision, parameter :: pi = acos(-1d0)

    entries(:) = 0
    rangeStarts(1) = 0
    do rangeIdx = 1, (nRanges-1)
        rangeStarts(rangeIdx+1) = rangeStarts(rangeIdx) + rangeReps(rangeIdx)
    enddo

    do rangeIdx = 1, nRanges
        range = ranges(rangeIdx)
        range2 = range**2
        !$omp parallel private(P, L1, L2, sweptAngle)
        !$omp do
        pointLoop: do pointIdx = (rangeStarts(rangeIdx) + 1), (rangeStarts(rangeIdx) + rangeReps(rangeIdx))
            sweptAngle = 0d0
            P = points(:,pointIdx)
            !check if the point is within range of any of the edges
            do vertIdx = 1, nVerts
                L1 = vertices(:,vertIdx)
                if (vertIdx == nVerts) then
                    L2 = vertices(:,1)
                else
                    L2 = vertices(:,vertIdx+1)
                endif
                if (MinSquaredDistance(L1-P, L2-P) <= range2) then
                    entries(pointIdx) = 1
                    cycle pointLoop
                endif
                sweptAngle = sweptAngle + SignedAngle(L1-P, L2-P)
            enddo
            if (abs(sweptAngle) > pi) then
                entries(pointIdx) = 2
            endif
        enddo pointLoop
        !$omp end do
        !$omp end parallel
    enddo
end subroutine PointNearPolygon