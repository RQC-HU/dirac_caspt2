! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

PROGRAM r4dcasci   ! DO CASCI CALC IN THIS PROGRAM!

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables
    use module_file_manager
    use module_2integrals
    use module_realonly
    use read_input_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer                    :: i0, nuniq, unit_eps, unit_input
    character(:), allocatable  :: filename

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

! MPI initialization and get the number of MPI processes (nprocs) and own process number.
#ifdef HAVE_MPI
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_rank(MPI_COMM_WORLD, rank, ierr)
#else
    rank = 0; nprocs = 1
#endif
    if (rank == 0) then
        print '(2(A,1X,I0))', 'initialization of mpi, rank :', rank, ' nprocs :', nprocs
        print *, ''
        print *, ' ENTER R4DCASCI PROGRAM written by M. Abe 2007.7.19'
        print *, ''
    end if
    tmem = 0.0d+00

    if (rank == 0) then
        call write_allocated_memory_size

        val = 0
        Call DATE_AND_TIME(VALUES=val)

        print *, 'Year = ', val(1), 'Mon = ', val(2), 'Date = ', val(3)
        print *, 'Hour = ', val(5), 'Min = ', val(6), 'Sec = ', val(7), '.', val(8)

        inittime = val(8)*(1.0d-03) + val(7) + val(6)*(6.0d+01) + val(5)*(6.0d+01)**2
        initdate = val(3)
        date0 = initdate; tsec0 = inittime

        print *, inittime
    end if

    call open_formatted_file(unit=unit_input, file='active.inp', status="old", optional_action='read')
    call read_input(unit_input)
    close (unit_input)

    if (rank == 0) then
        print *, 'ninact        =', ninact
        print *, 'nact          =', nact
        print *, 'nsec          =', nsec
        print *, 'nelec         =', nelec
        print *, 'nroot         =', nroot
        print *, 'selectroot    =', selectroot
        print *, 'totsym        =', totsym
        print *, 'eshift        =', eshift
        print *, 'dirac_version =', dirac_version
        if (ras1_size /= 0) print *, "RAS1 =", ras1_list
        if (ras2_size /= 0) print *, "RAS2 =", ras2_list
        if (ras3_size /= 0) print *, "RAS3 =", ras3_list
    end if

    ! Read MRCONEE file (orbital energies, symmetries and multiplication tables)
    filename = 'MRCONEE'
    call read_mrconee(filename)

    ! Read around the MDCINT file and determine if the imaginary part of the 2-electron integral is written or not.
    call check_realonly()
    if (skip_mdcint) then
        if (rank == 0) print *, "Skip create_newmdcint (Activated skip_mdcint option by user input file)"
    else
        call timing(date0, tsec0, date1, tsec1)
        date0 = date1; tsec0 = tsec1
        ! Create UTChem type MDCINT file from Dirac MDCINT file
        call create_newmdcint
        call timing(date0, tsec0, date1, tsec1)
        date0 = date1; tsec0 = tsec1
    end if

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

#ifdef DEBUG
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
#endif

!! NOW MAKE FOCK MATRIX FOR CASCI STATE
!! fij = hij + SIGUMA_kl[<0|Ekl|0>{(ij|kl)-(il|kj)}
    if (rank == 0) then
        call timing(date0, tsec0, date1, tsec1)
        date0 = date1; tsec0 = tsec1
    end if
    if (realonly%is_realonly()) then
        if (.not. allocated(fock_real)) then
            allocate (fock_real(nmo, nmo)); call memplus(KIND(fock_real), SIZE(fock_real), 1)
        end if
        fock_real(:, :) = 0.0d+00
        Call fockcasci_real
    else
        if (.not. allocated(fock_cmplx)) then
            allocate (fock_cmplx(nmo, nmo)); call memplus(KIND(fock_cmplx), SIZE(fock_cmplx), 2)
        end if
        fock_cmplx(:, :) = 0.0d+00
        Call fockcasci_complex
    end if

    if (rank == 0) then
        print *, 'end building fock'
        call timing(date0, tsec0, date1, tsec1)
        date0 = date1; tsec0 = tsec1
    end if

#ifdef DEBUG
    call prtoutfock
#endif

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

    ! Deallocate the memory
    if (allocated(ras1_list)) then
        Call memminus(KIND(ras1_list), SIZE(ras1_list), 1); deallocate (ras1_list)
    end if
    if (allocated(ras2_list)) then
        Call memminus(KIND(ras2_list), SIZE(ras2_list), 1); deallocate (ras2_list)
    end if
    if (allocated(ras3_list)) then
        Call memminus(KIND(ras3_list), SIZE(ras3_list), 1); deallocate (ras3_list)
    end if
    if (allocated(space_idx)) then
        Call memminus(KIND(space_idx), SIZE(space_idx), 1); deallocate (space_idx)
    end if
    if (allocated(cir)) then
        Call memminus(KIND(cir), SIZE(cir), 1); deallocate (cir)
    end if
    if (allocated(cii)) then
        Call memminus(KIND(cii), SIZE(cii), 1); deallocate (cii)
    end if
    if (allocated(eigen)) then
        Call memminus(KIND(eigen), SIZE(eigen), 1); deallocate (eigen)
    end if
    if (allocated(fock_real)) then
        Call memminus(KIND(fock_real), SIZE(fock_real), 1); deallocate (fock_real)
    end if
    if (allocated(fock_cmplx)) then
        Call memminus(KIND(fock_cmplx), SIZE(fock_cmplx), 2); deallocate (fock_cmplx)
    end if
    if (allocated(eps)) then
        Call memminus(KIND(eps), SIZE(eps), 1); deallocate (eps)
    end if
    if (allocated(MULTB_S)) then
        Call memminus(KIND(MULTB_S), SIZE(MULTB_S), 1); deallocate (MULTB_S)
    end if
    if (allocated(MULTB_D)) then
        Call memminus(KIND(MULTB_D), SIZE(MULTB_D), 1); deallocate (MULTB_D)
    end if
    if (allocated(MULTB_DS)) then
        Call memminus(KIND(MULTB_DS), SIZE(MULTB_DS), 1); deallocate (MULTB_DS)
    end if
    if (allocated(irpamo)) then
        Call memminus(KIND(irpamo), SIZE(irpamo), 1); deallocate (irpamo)
    end if
    if (allocated(indmo_cas_to_dirac)) then
        Call memminus(KIND(indmo_cas_to_dirac), SIZE(indmo_cas_to_dirac), 1); deallocate (indmo_cas_to_dirac)
    end if
    if (allocated(indmo_dirac_to_cas)) then
        Call memminus(KIND(indmo_dirac_to_cas), SIZE(indmo_dirac_to_cas), 1); deallocate (indmo_dirac_to_cas)
    end if
    if (allocated(one_elec_int_i)) then
        Call memminus(KIND(one_elec_int_i), SIZE(one_elec_int_i), 1); deallocate (one_elec_int_i)
    end if
    if (allocated(inttwi)) then
        Call memminus(KIND(inttwi), SIZE(inttwi), 1); deallocate (inttwi)
    end if
    if (allocated(one_elec_int_r)) then
        Call memminus(KIND(one_elec_int_r), SIZE(one_elec_int_r), 1); deallocate (one_elec_int_r)
    end if
    if (allocated(inttwr)) then
        Call memminus(KIND(inttwr), SIZE(inttwr), 1); deallocate (inttwr)
    end if
    if (allocated(int2r_f1)) then
        Call memminus(KIND(int2r_f1), SIZE(int2r_f1), 1); deallocate (int2r_f1)
    end if
    if (allocated(int2i_f1)) then
        Call memminus(KIND(int2i_f1), SIZE(int2i_f1), 1); deallocate (int2i_f1)
    end if
    if (allocated(int2r_f2)) then
        Call memminus(KIND(int2r_f2), SIZE(int2r_f2), 1); deallocate (int2r_f2)
    end if
    if (allocated(int2i_f2)) then
        Call memminus(KIND(int2i_f2), SIZE(int2i_f2), 1); deallocate (int2i_f2)
    end if
    if (rank == 0) then
        call write_allocated_memory_size

        Call timing(val(3), totalsec, date0, tsec0)
        print *, 'End r4dcasci part'
    end if
#ifdef HAVE_MPI
    call MPI_FINALIZE(ierr)
#endif

END program r4dcasci
