! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE prtoutfock  ! TO PRINT OUT FOCK MATRIX

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables

    Implicit NONE
    integer :: i, j

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    print *, 'inactive-inactive'

    do i = 1, ninact
        do j = 1, ninact
            if ((i /= j) .and. (ABS(fock_cmplx(i, j)) > 1.0d-10)) then
                print '(2I4,3E20.10)', i, j, fock_cmplx(i, j), one_elec_int_r(i, j)
            end if
        end do
    end do

    print *, 'inactive-active'

    do i = 1, ninact
        do j = ninact + 1, ninact + nact
            if ((i /= j) .and. (ABS(fock_cmplx(i, j)) > 1.0d-10)) then
                print '(2I4,3E20.10)', i, j, fock_cmplx(i, j), one_elec_int_r(i, j)
            end if
        end do
    end do

    print *, 'inactive-secondary'

    do i = 1, ninact
        do j = ninact + nact + 1, ninact + nact + nsec
            if ((i /= j) .and. (ABS(fock_cmplx(i, j)) > 1.0d-10)) then
                print '(2I4,3E20.10)', i, j, fock_cmplx(i, j), one_elec_int_r(i, j)
            end if
        end do
    end do

    print *, 'active-active'

    do i = ninact + 1, ninact + nact
        do j = ninact + 1, ninact + nact
            if ((i /= j) .and. (ABS(fock_cmplx(i, j)) > 1.0d-10)) then
                print '(2I4,3E20.10)', i, j, fock_cmplx(i, j), one_elec_int_r(i, j)
            end if
        end do
    end do

    print *, 'active-secondary'

    do i = ninact + 1, ninact + nact
        do j = ninact + nact + 1, ninact + nact + nsec
            if ((i /= j) .and. (ABS(fock_cmplx(i, j)) > 1.0d-10)) then
                print '(2I4,3E20.10)', i, j, fock_cmplx(i, j), one_elec_int_r(i, j)
            end if
        end do
    end do

    print *, 'secondary-secondary'

    do i = ninact + nact + 1, ninact + nact + nsec
        do j = ninact + nact + 1, ninact + nact + nsec
            if ((i /= j) .and. (ABS(fock_cmplx(i, j)) > 1.0d-10)) then
                print '(2I4,3E20.10)', i, j, fock_cmplx(i, j), one_elec_int_r(i, j)
            end if
        end do
    end do
end SUBROUTINE prtoutfock
