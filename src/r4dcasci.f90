! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

PROGRAM r4dcasci   ! DO CASCI CALC IN THIS PROGRAM!

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module
    use module_file_manager
    use module_2integrals
    use module_realonly
    use read_input_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer                 :: i0, nuniq, inisym, endsym, unit_eps, unit_input
    character*50            :: filename

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
    debug = .FALSE.
    tmem = 0.0d+00

    if (rank == 0) then
        print '("Current Memory is ",F10.2,"MB")', tmem/1024/1024

        val = 0
        Call DATE_AND_TIME(VALUES=val)

        print *, 'Year = ', val(1), 'Mon = ', val(2), 'Date = ', val(3)
        print *, 'Hour = ', val(5), 'Min = ', val(6), 'Sec = ', val(7), '.', val(8)

        totalsec = val(8)*(1.0d-03) + val(7) + val(6)*(6.0d+01) + val(5)*(6.0d+01)**2
        initdate = val(3)
        inittime = totalsec

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
        print *, 'ncore         =', ncore
        print *, 'nbas          =', nbas
        print *, 'eshift        =', eshift
        print *, 'dirac_version =', dirac_version
        if (ras1_size /= 0) print *, "RAS1 =", ras1_list
        if (ras2_size /= 0) print *, "RAS2 =", ras2_list
        if (ras3_size /= 0) print *, "RAS3 =", ras3_list
    end if

    ! Read MRCONEE file (orbital energies, symmetries and multiplication tables)
    filename = 'MRCONEE'
    call readorb_enesym(filename)
    call read1mo(filename)

    if (rank == 0) print *, 'realc', realc, ECORE, ninact, nact, nsec, nmo

    call check_realonly()
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
    if (rank == 0) print *, "iwamuro modify"
    If (mod(nelec, 2) == 0) then
        inisym = nsymrpa + 1
        endsym = 2*nsymrpa
    Else
        inisym = 1
        endsym = nsymrpa
    End if

    ! Print the irreducible representation used to calculate CASCI energy.
    if (rank == 0) then
        print '("Current Memory is ",F10.2,"MB")', tmem/1024/1024

        print *, ' '
        print *, '*******************************'
        print *, ' '
        print *, 'IREP IS ', repna(totsym)
        print *, ' '
        print *, '*******************************'
        print *, ' '
    end if
    realcvec = .TRUE.

    ! Calculate 0th order energy (CASCI energy)
    Call casci

    realc = .FALSE.      !!!      realc =.TRUE.
    realcvec = .FALSE.   !!!      realcvec =.TRUE.

    if (rank == 0) then
        print *, 'FOR TEST WE DO (F,F)'
        print *, realc, 'realc'
        print *, realcvec, 'realcvec'
    end if
    iroot = selectroot

    Call e0test

    if (rank == 0) print '("Current Memory is ",F10.2,"MB")', tmem/1024/1024
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
        call fockhf1_real
    else
        Allocate (fock_cmplx(nmo, nmo)); Call memplus(KIND(fock_cmplx), SIZE(fock_cmplx), 2)
        fock_cmplx(:, :) = 0.0d+00
        call fockhf1
    End if
    ! debug = .TRUE.
    ! If (debug) then
    !     if (rank == 0) print *, 'fockhf1 start'
    !     Call fockhf1
    ! End if

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
        Call fockcasci
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

    ! Deallocate the memory
    if (allocated(ras1_list)) deallocate (ras1_list); Call memminus(KIND(ras1_list), SIZE(ras1_list), 1)
    if (allocated(ras2_list)) deallocate (ras2_list); Call memminus(KIND(ras2_list), SIZE(ras2_list), 1)
    if (allocated(ras3_list)) deallocate (ras3_list); Call memminus(KIND(ras3_list), SIZE(ras3_list), 1)
    if (allocated(space_idx)) deallocate (space_idx); Call memminus(KIND(space_idx), SIZE(space_idx), 1)
    if (allocated(cir)) deallocate (cir); Call memminus(KIND(cir), SIZE(cir), 1)
    if (allocated(cii)) deallocate (cii); Call memminus(KIND(cii), SIZE(cii), 1)
    if (allocated(eigen)) deallocate (eigen); Call memminus(KIND(eigen), SIZE(eigen), 1)
    if (allocated(fock_real)) deallocate (fock_real); Call memminus(KIND(fock_real), SIZE(fock_real), 1)
    if (allocated(fock_cmplx)) deallocate (fock_cmplx); Call memminus(KIND(fock_cmplx), SIZE(fock_cmplx), 2)
    if (allocated(eps)) deallocate (eps); Call memminus(KIND(eps), SIZE(eps), 1)
    if (allocated(cas_idx)) deallocate (cas_idx); Call memminus(KIND(cas_idx), SIZE(cas_idx), 1)
    if (allocated(MULTB_S)) deallocate (MULTB_S); Call memminus(KIND(MULTB_S), SIZE(MULTB_S), 1)
    if (allocated(MULTB_D)) deallocate (MULTB_D); Call memminus(KIND(MULTB_D), SIZE(MULTB_D), 1)
    if (allocated(MULTB_DS)) deallocate (MULTB_DS); Call memminus(KIND(MULTB_DS), SIZE(MULTB_DS), 1)
    if (allocated(MULTB_DF)) deallocate (MULTB_DF); Call memminus(KIND(MULTB_DF), SIZE(MULTB_DF), 1)
    if (allocated(MULTB_DB)) deallocate (MULTB_DB); Call memminus(KIND(MULTB_DB), SIZE(MULTB_DB), 1)
    if (allocated(MULTB_SB)) deallocate (MULTB_SB); Call memminus(KIND(MULTB_SB), SIZE(MULTB_SB), 1)
    if (allocated(irpamo)) deallocate (irpamo); Call memminus(KIND(irpamo), SIZE(irpamo), 1)
    if (allocated(indmo_cas_to_dirac)) then
        deallocate (indmo_cas_to_dirac); Call memminus(KIND(indmo_cas_to_dirac), SIZE(indmo_cas_to_dirac), 1)
    end if
    if (allocated(indmo_dirac_to_cas)) then
        deallocate (indmo_dirac_to_cas); Call memminus(KIND(indmo_dirac_to_cas), SIZE(indmo_dirac_to_cas), 1)
    end if
    if (allocated(one_elec_int_i)) deallocate (one_elec_int_i); Call memminus(KIND(one_elec_int_i), SIZE(one_elec_int_i), 1)
    if (allocated(inttwi)) deallocate (inttwi); Call memminus(KIND(inttwi), SIZE(inttwi), 1)
    if (allocated(one_elec_int_r)) deallocate (one_elec_int_r); Call memminus(KIND(one_elec_int_r), SIZE(one_elec_int_r), 1)
    if (allocated(inttwr)) deallocate (inttwr); Call memminus(KIND(inttwr), SIZE(inttwr), 1)
    if (allocated(int2r_f1)) deallocate (int2r_f1); Call memminus(KIND(int2r_f1), SIZE(int2r_f1), 1)
    if (allocated(int2i_f1)) deallocate (int2i_f1); Call memminus(KIND(int2i_f1), SIZE(int2i_f1), 1)
    if (allocated(int2r_f2)) deallocate (int2r_f2); Call memminus(KIND(int2r_f2), SIZE(int2r_f2), 1)
    if (allocated(int2i_f2)) deallocate (int2i_f2); Call memminus(KIND(int2i_f2), SIZE(int2i_f2), 1)
    if (rank == 0) then
        print '("Current Memory is ",F10.2,"MB")', tmem/1024/1024

        Call timing(val(3), totalsec, date0, tsec0)
        print *, 'End r4dcasci part'
    end if
#ifdef HAVE_MPI
    call MPI_FINALIZE(ierr)
#endif

END program r4dcasci