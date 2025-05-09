! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE fockdiag

! Diagonalize the Fock matrix.

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_file_manager
    use module_realonly, only: realonly
    use module_global_variables

    Implicit NONE

    integer                 :: i0, n, n0, n1, nspace(3, 3)
    real(8), allocatable :: fa(:, :)
    complex*16, allocatable :: fac(:, :)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    if (debug .and. rank == 0) print *, 'fockdiag start'

    If (realonly%is_realonly()) then          ! real(8)
        Allocate (fa(nmo, nmo)); Call memplus(KIND(fa), SIZE(fa), 1)
        eps = 0.0d+00
        fa = 0.0d+00
    Else
        Allocate (fac(nmo, nmo)); Call memplus(KIND(fac), SIZE(fac), 2)
        eps = 0.0d+00
        fac = 0.0d+00
    End if

    nspace(1, 1) = 1
    nspace(2, 1) = ninact
    nspace(3, 1) = ninact

    nspace(1, 2) = global_act_start
    nspace(2, 2) = global_act_end
    nspace(3, 2) = nact

    nspace(1, 3) = global_sec_start
    nspace(2, 3) = global_sec_end
    nspace(3, 3) = nsec

    Do i0 = 1, 3
        n0 = nspace(1, i0)
        n1 = nspace(2, i0)
        n = nspace(3, i0)
        if (debug .and. rank == 0) then
            if (i0 == 1) print *, 'FOR INACTIVE-INACTIVE ROTATION !'
            if (i0 == 2) print *, 'FOR ACTIVE-ACTIVE ROTATION !'
            if (i0 == 3) print *, 'FOR SECONDARY-SECONDARY ROTATION !'
        end if
        if (n == 0) cycle
        if (realonly%is_realonly()) then
            Call rdiag0(n, n0, n1, fa(n0:n1, n0:n1), eps(n0:n1))
        else
            Call cdiag0(n, n0, n1, fac(n0:n1, n0:n1), eps(n0:n1))
        end if
    End do

    if (realonly%is_realonly()) then
        Call traci(fa(global_act_start:global_act_end, global_act_start:global_act_end))
        fock_real(1:nmo, 1:nmo) = fa(1:nmo, 1:nmo)
        Call memminus(KIND(fa), SIZE(fa), 1); deallocate (fa)
    else
        Call tracic(fac(global_act_start:global_act_end, global_act_start:global_act_end))
        fock_cmplx(1:nmo, 1:nmo) = fac(1:nmo, 1:nmo)
        Call memminus(KIND(fac), SIZE(fac), 2); deallocate (fac)
    end if
#ifdef DEBUG
    Call e0aftertra
#endif
    if (debug .and. rank == 0) print *, 'fockdiag end'
end subroutine fockdiag
