! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

subroutine r4dcasci   ! DO CASCI CALC IN THIS PROGRAM!

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables
    ! use module_file_manager, only: open_unformatted_file
    use module_file_manager
    use module_2integrals
    use module_realonly
    use read_input_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer                    :: i0, nuniq, unit_eps

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    if (rank == 0) then
        call write_allocated_memory_size

        val = 0
        Call DATE_AND_TIME(VALUES=val)

        print *, 'Year = ', val(1), 'Mon = ', val(2), 'Date = ', val(3)
        print *, 'Hour = ', val(5), 'Min = ', val(6), 'Sec = ', val(7), '.', val(8)

        totalsec = val(8)*(1.0d-03) + val(7) + val(6)*(6.0d+01) + val(5)*(6.0d+01)**2
        initdate = val(3)
        inittime = totalsec

        print *, inittime
    end if

    if (skip_mdcint) then
        if (rank == 0) print *, "Skip create_newmdcint (Activated skip_mdcint option by user input file)"
    else
        ! Create UTChem type MDCINT file from Dirac MDCINT file
        call create_newmdcint
    end if
    if (rank == 0) print '(a)', 'Before readint2_casci'

    ! Read UTChem type MDCINT files and expands the 2-electron integral in memory
    Call readint2_casci(mdcintnew, nuniq)

    if (rank == 0) print *, 'nmo        =', nmo
    nmo = ninact + nact + nsec

    ! Print the irreducible representation used to calculate CASCI energy.
    if (rank == 0) then
        call write_allocated_memory_size
        print *, ' '
        print *, '*******************************'
        print *, ' '
        print *, 'IREP IS ', repna(totsym)
        print *, ' '
        print *, '*******************************'
        print *, ' '
    end if

    ! Calculate 0th order energy (CASCI energy)
    Call casci
    iroot = selectroot

    ! Recalculate the 0th order energy (CASCI energy) using the 1,2 electron integrals adn CI coefficients
    Call e0test

    call write_allocated_memory_size
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!            BUILDING  FOCK MATRIX               !
!  fij = hij + SIGUMA[<0|Ekl|0>{(ij|kl)-(il|kj)} !
!                 kl                             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

!! TEST TO CALCULATE FOCK MATRIX OF HF STATE fpq = hpq + SIGUMA_r[(pq|rr)-(pr|qr)]
!! THIS MUST BE DIAGONAL MATRIX AND DIAGONAL ELEMENTS CORESPONDS TO SPINOR ENERGIES.
    if (realonly%is_realonly()) then
        allocate (fock_real(nmo, nmo)); Call memplus(KIND(fock_real), SIZE(fock_real), 1)
        fock_real(:, :) = 0.0d+00
        call fock_matrix_of_hf_real
    else
        Allocate (fock_cmplx(nmo, nmo)); Call memplus(KIND(fock_cmplx), SIZE(fock_cmplx), 2)
        fock_cmplx(:, :) = 0.0d+00
        call fock_matrix_of_hf_complex
    End if

!! NOW MAKE FOCK MATRIX FOR CASCI STATE
!! fij = hij + SIGUMA_kl[<0|Ekl|0>{(ij|kl)-(il|kj)}
    if (rank == 0) then
        print *, 'before building fock'
        date1 = date0
        tsec1 = tsec0
        Call timing(date1, tsec1, date0, tsec0)
    end if
    if (realonly%is_realonly()) then
        fock_real(:, :) = 0.0d+00
        Call fockcasci_real
    else
        fock_cmplx(:, :) = 0.0d+00
        Call fockcasci_complex
    end if

    if (rank == 0) then
        print *, 'end building fock'
        date1 = date0
        tsec1 = tsec0
        Call timing(date1, tsec1, date0, tsec0)
    end if
    debug = .FALSE.
    if (rank == 0) print *, debug, 'debug'
    if (debug) Call prtoutfock

    Allocate (eps(nmo)); Call memplus(KIND(eps), SIZE(eps), 1)
    eps = 0.0d+00

    ! Diagonalize the Fock matrix
    Call fockdiag

    ! Print orbital energies
    if (rank == 0) then
        Do i0 = 1, nmo
            print *, 'eps(', i0, ')=', eps(i0)
        End do
    end if

    ! Store orbital energies in EPS file
    if (rank == 0) then ! Only master ranks are allowed to create files used by CASPT2 except for MDCINTNEW.
        call open_unformatted_file(unit=unit_eps, file="EPS", status="replace", optional_action="write")
        write (unit_eps) nmo
        write (unit_eps) eps(1:nmo)
        close (unit_eps)
    end if

    if (rank == 0) then
        call write_allocated_memory_size

        Call timing(val(3), totalsec, date0, tsec0)
        print *, 'End r4dcasci part'
    end if

end subroutine r4dcasci
