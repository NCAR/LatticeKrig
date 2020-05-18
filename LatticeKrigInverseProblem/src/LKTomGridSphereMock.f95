subroutine LKTomGridSphereMock(vals, nVals)
    implicit none
    integer, intent(in) :: nVals
    integer :: idx

    double precision :: vals(nVals), compare(nVals)

    compare(:) = 2
    vals(1) = sum(vals * compare)
!    do idx = 1, nVals
!        vals(idx) = vals(idx) * vals(idx)
!    enddo
end subroutine LKTomGridSphereMock