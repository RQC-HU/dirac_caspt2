! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE prtoutfock  ! TO PRINT OUT FOCK MATRIX

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables
    use module_realonly, only: realonly

    Implicit NONE
    integer :: i, j
    real(8), parameter :: prtfock_threshold = 1.0e-10

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    if (realonly%is_realonly()) then
        call prtoutfock_realonly
    else
        call prtoutfock_cmplx
    end if

contains

    subroutine prtoutfock_cmplx
        implicit none
        print *, 'inactive-inactive'

        do i = 1, ninact
            do j = 1, ninact
                if ((i /= j) .and. (ABS(fock_cmplx(i, j)) > prtfock_threshold)) then
                    print '(2I4,3E20.10)', i, j, fock_cmplx(i, j), one_elec_int_r(i, j)
                end if
            end do
        end do

        print *, 'inactive-active'

        do i = 1, ninact
            do j = global_act_start, global_act_end
                if ((i /= j) .and. (ABS(fock_cmplx(i, j)) > prtfock_threshold)) then
                    print '(2I4,3E20.10)', i, j, fock_cmplx(i, j), one_elec_int_r(i, j)
                end if
            end do
        end do

        print *, 'inactive-secondary'

        do i = 1, ninact
            do j = global_sec_start, global_sec_end
                if ((i /= j) .and. (ABS(fock_cmplx(i, j)) > prtfock_threshold)) then
                    print '(2I4,3E20.10)', i, j, fock_cmplx(i, j), one_elec_int_r(i, j)
                end if
            end do
        end do

        print *, 'active-active'

        do i = global_act_start, global_act_end
            do j = global_act_start, global_act_end
                if ((i /= j) .and. (ABS(fock_cmplx(i, j)) > prtfock_threshold)) then
                    print '(2I4,3E20.10)', i, j, fock_cmplx(i, j), one_elec_int_r(i, j)
                end if
            end do
        end do

        print *, 'active-secondary'

        do i = global_act_start, global_act_end
            do j = global_sec_start, global_sec_end
                if ((i /= j) .and. (ABS(fock_cmplx(i, j)) > prtfock_threshold)) then
                    print '(2I4,3E20.10)', i, j, fock_cmplx(i, j), one_elec_int_r(i, j)
                end if
            end do
        end do

        print *, 'secondary-secondary'

        do i = global_sec_start, global_sec_end
            do j = global_sec_start, global_sec_end
                if ((i /= j) .and. (ABS(fock_cmplx(i, j)) > prtfock_threshold)) then
                    print '(2I4,3E20.10)', i, j, fock_cmplx(i, j), one_elec_int_r(i, j)
                end if
            end do
        end do
    end subroutine prtoutfock_cmplx

    subroutine prtoutfock_realonly
        implicit none
        print *, 'inactive-inactive'

        do i = 1, ninact
            do j = 1, ninact
                if ((i /= j) .and. (ABS(fock_real(i, j)) > prtfock_threshold)) then
                    print '(2I4,2E20.10)', i, j, fock_real(i, j), one_elec_int_r(i, j)
                end if
            end do
        end do

        print *, 'inactive-active'

        do i = 1, ninact
            do j = global_act_start, global_act_end
                if ((i /= j) .and. (ABS(fock_real(i, j)) > prtfock_threshold)) then
                    print '(2I4,2E20.10)', i, j, fock_real(i, j), one_elec_int_r(i, j)
                end if
            end do
        end do

        print *, 'inactive-secondary'

        do i = 1, ninact
            do j = global_sec_start, global_sec_end
                if ((i /= j) .and. (ABS(fock_real(i, j)) > prtfock_threshold)) then
                    print '(2I4,2E20.10)', i, j, fock_real(i, j), one_elec_int_r(i, j)
                end if
            end do
        end do

        print *, 'active-active'

        do i = global_act_start, global_act_end
            do j = global_act_start, global_act_end
                if ((i /= j) .and. (ABS(fock_real(i, j)) > prtfock_threshold)) then
                    print '(2I4,2E20.10)', i, j, fock_real(i, j), one_elec_int_r(i, j)
                end if
            end do
        end do

        print *, 'active-secondary'

        do i = global_act_start, global_act_end
            do j = global_sec_start, global_sec_end
                if ((i /= j) .and. (ABS(fock_real(i, j)) > prtfock_threshold)) then
                    print '(2I4,2E20.10)', i, j, fock_real(i, j), one_elec_int_r(i, j)
                end if
            end do
        end do

        print *, 'secondary-secondary'

        do i = global_sec_start, global_sec_end
            do j = global_sec_start, global_sec_end
                if ((i /= j) .and. (ABS(fock_real(i, j)) > prtfock_threshold)) then
                    print '(2I4,2E20.10)', i, j, fock_real(i, j), one_elec_int_r(i, j)
                end if
            end do
        end do
    end subroutine prtoutfock_realonly

end SUBROUTINE prtoutfock
