
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE traci(fa)  ! Transform CI matrix for new spinor basis

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_error, only: stop_with_errorcode
    use module_global_variables
    use module_file_manager, only: open_unformatted_file

    Implicit NONE
    real(8), intent(in)  :: fa(global_act_start:global_act_end, global_act_start:global_act_end)

    integer :: i0, j0, i, info
    integer :: ok, unit_newcicoeff
    integer :: occ(nelec, ndet)

    integer, allocatable    :: IPIV(:)
    real(8), allocatable    :: ds(:, :)
    complex*16, allocatable ::  ci(:)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    occ = 0
    if (rank == 0) print *, 'Enter TRACI'

    Do i0 = 1, ndet
        i = 0
        ok = 0
        Do j0 = 0, 63 ! 64 bits integer are possible with 64 spinors
            if (btest(get_val(dict_cas_idx, i0), j0)) then ! This condition should be true nelec times
                i = i + 1
                if (j0 + 1 <= nact) then ! j0+1 means occupied spinor labeled by casci
                    occ(i, i0) = j0 + 1  ! This is energetic order inside active spinor!
                    ok = ok + 1
                End if
            end if
        End do
    End do

    Allocate (ds(ndet, ndet))
    ! ds = <k|k~>
    Do i0 = 1, ndet     ! k  (old)
        Do j0 = 1, ndet  ! k~ (new)   <k|k~>
            Call dets(fa(global_act_start:global_act_end, global_act_start:global_act_end), &
                      occ(1:nelec, i0), occ(1:nelec, j0), ds(i0, j0))
        End do
    End do

    Allocate (IPIV(ndet))
    ! We want to create a cir vector rotated by ds^(-1) matrix
    ! rotated_ci = ds^(-1) * cir
    ! Therefore we need to solve ds * rotated_ci = cir, so we can get rotated_ci by simply calling dgesv.
    ! dgesv: https://netlib.org/lapack/explore-html/d7/d3b/group__double_g_esolve_ga5ee879032a8365897c3ba91e3dc8d512.html
    ! ds is overwritten by its LU decomposition
    ! cir is overwritten by rotated_ci
    if (rank == 0) print *, 'Solve linear equations for ci matrix using dgesv'
    Call dgesv(ndet, 1, ds, ndet, IPIV, cir(1:ndet, selectroot), ndet, info)
    Deallocate (IPIV)
    if (info /= 0) then
        if (rank == 0) print *, 'Error in dgesv, info = ', info
        call stop_with_errorcode(info)
    end if

    allocate (ci(ndet)); call memplus(kind(ci), size(ci), 2)
    ci(:) = (0.0d+00, 0.0d+00)
    ci = cir(1:ndet, selectroot) ! Imaginary part of ci is zero
    if (rank == 0) then ! Only master ranks are allowed to create files used by CASPT2 except for MDCINTNEW.
        call open_unformatted_file(unit=unit_newcicoeff, file="NEWCICOEFF", status='replace', optional_action='write')
        write (unit_newcicoeff) ci(:)
        close (unit_newcicoeff)
    end if
    if (allocated(ci)) call memminus(kind(ci), size(ci), 2); Deallocate (ci)

    Deallocate (ds)

End subroutine traci

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE tracic(fac)  ! Transform CI matrix for new spinor basis

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_error, only: stop_with_errorcode
    use module_global_variables
    use module_file_manager, only: open_unformatted_file
    use module_time

    Implicit NONE

    complex*16, intent(in)  :: fac(global_act_start:global_act_end, global_act_start:global_act_end)

    integer :: i0, j0, i, info
    integer :: ok, unit_newcicoeff
    integer :: occ(nelec, ndet)

    integer, allocatable     :: IPIV(:)
    complex*16, Allocatable  :: ds(:, :), ci(:)
    type(time_type)          :: tmp_start_time, tmp_end_time
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    occ = 0
    if (rank == 0) print *, 'Enter TRACI'
    call get_current_time(tmp_start_time)

    Do i0 = 1, ndet
        i = 0
        ok = 0
        Do j0 = 0, 63 ! 64 bits integer are possible with 64 spinors
            if (btest(get_val(dict_cas_idx, i0), j0)) then ! This condition should be true nelec times
                i = i + 1
                if (j0 + 1 <= nact) then ! j0+1 means occupied spinor labeled by casci
                    occ(i, i0) = j0 + 1  ! This is energetic order inside active spinor!
                    ok = ok + 1
                End if
            end if
        End do
    End do

    Allocate (ds(ndet, ndet))
    Do i0 = 1, ndet     ! k  (old)
        Do j0 = 1, ndet  ! k~ (new)   <k|k~>
            Call detsc(fac(global_act_start:global_act_end, global_act_start:global_act_end), &
                       occ(1:nelec, i0), occ(1:nelec, j0), ds(i0, j0))
        End do
    End do
    if (rank == 0) print *, 'End detsc'

    ! Preparation for solving linear equations for ci matrix using zgesv
    call get_current_time_and_print_diff(tmp_start_time, tmp_end_time); tmp_start_time = tmp_end_time
    Allocate (IPIV(ndet))
    Allocate (ci(ndet))
    ci = DCMPLX(cir(1:ndet, selectroot), cii(1:ndet, selectroot))

    ! We want to create a ci vector rotated by ds^(-1) matrix
    ! rotated_ci = ds^(-1) * ci
    ! Therefore we need to solve ds * rotated_ci = ci, so we can get rotated_ci by simply calling zgesv.
    ! zgesv: https://netlib.org/lapack/explore-html/d6/d10/group__complex16_g_esolve_ga531713dfc62bc5df387b7bb486a9deeb.html
    ! ds is overwritten by its LU decomposition
    ! ci is overwritten by rotated_ci
    if (rank == 0) print *, 'Solve linear equations for ci matrix using zgesv'
    Call zgesv(ndet, 1, ds, ndet, IPIV, ci, ndet, info) ! ds is overwritten by its LU decomposition, ci is overwritten by rotated_ci
    Deallocate (IPIV)
    if (info /= 0) then
        if (rank == 0) print *, 'Error in zgesv, info = ', info
        call stop_with_errorcode(info)
    end if
    if (rank == 0) print *, 'End zgesv', rank
    call get_current_time_and_print_diff(tmp_start_time, tmp_end_time)
    ! ci is now rotated_ci and we need to rotate it back to cir and cii
    cir(1:ndet, selectroot) = DBLE(ci(1:ndet))
    cii(1:ndet, selectroot) = DIMAG(ci(1:ndet))
    if (rank == 0) then ! Only master ranks are allowed to create files used by CASPT2 except for MDCINTNEW.
        call open_unformatted_file(unit=unit_newcicoeff, file="NEWCICOEFF", status='replace', optional_action='write')
        write (unit_newcicoeff) ci(1:ndet)
        close (unit_newcicoeff)
    end if

    Deallocate (ci)
    Deallocate (ds)

End subroutine tracic
