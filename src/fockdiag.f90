! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE fockdiag

! Diagonalize the Fock matrix.

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_file_manager
    use module_realonly, only: realonly
    use four_caspt2_module

    Implicit NONE

    integer                 ::  i, j
    integer                 :: unit_transfock
    integer                 :: i0, n, n0, n1, nspace(3, 3)
    real*8, allocatable :: fa(:, :)
    complex*16, allocatable :: fac(:, :)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    if (rank == 0) print *, 'fockdiag start'

    If (realonly%is_realonly()) then          ! real*8
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

    nspace(1, 2) = ninact + 1
    nspace(2, 2) = ninact + nact
    nspace(3, 2) = nact

    nspace(1, 3) = ninact + nact + 1
    nspace(2, 3) = ninact + nact + nsec
    nspace(3, 3) = nsec

    Do i0 = 1, 3

        n0 = nspace(1, i0)
        n1 = nspace(2, i0)
        n = nspace(3, i0)

        if (rank == 0) then
            if (i0 == 1) print *, 'FOR INACTIVE-INACTIVE ROTATION !'
            if (i0 == 2) print *, 'FOR ACTIVE-ACTIVE ROTATION !'
            if (i0 == 3) print *, 'FOR SECONDARY-SECONDARY ROTATION !'
        end if
        if (realonly%is_realonly()) then

            Call rdiag0(n, n0, n1, fa(n0:n1, n0:n1), eps(n0:n1))

            ! write (5) n0, n1, n
            ! write (5) fa(n0:n1, n0:n1)
            if (rank == 0) then
                print *, n0, n1, n

                print *, 'fa '

                do i = n0, n1
                    print '(30E13.5)', (fa(i, j), j=n0, n1)
                end do

                print *, 'fock_real'
                do i = n0, n1
                    print '(30E13.5)', (fock_real(i, j), j=n0, n1)
                end do
            end if
        else
            Call cdiag0(n, n0, n1, fac(n0:n1, n0:n1), eps(n0:n1))

        end if

    End do ! i0

    if (realonly%is_realonly()) then

        Call traci(fa(ninact + 1:ninact + nact, ninact + 1:ninact + nact))

        fock_real(1:nmo, 1:nmo) = fa(1:nmo, 1:nmo)

        Call e0aftertra

        deallocate (fa); Call memminus(KIND(fa), SIZE(fa), 1)

    else

        Call tracic(fac(ninact + 1:ninact + nact, ninact + 1:ninact + nact))

        fock_cmplx(1:nmo, 1:nmo) = fac(1:nmo, 1:nmo)

        Call e0aftertrac

        deallocate (fac); Call memminus(KIND(fac), SIZE(fac), 2)

    end if

    if (rank == 0) then ! Only master ranks are allowed to create files used by CASPT2 except for MDCINTNEW.
        call open_unformatted_file(unit=unit_transfock, file='TRANSFOCK', status='replace', optional_action='write')
        write (unit_transfock) nmo
        if (realonly%is_realonly()) then
            write (unit_transfock) fock_real(1:nmo, 1:nmo)
        else
            write (unit_transfock) fock_cmplx(1:nmo, 1:nmo)
        end if
        close (unit_transfock)
    end if

    if (rank == 0) print *, 'fockdiag end'
end subroutine fockdiag
