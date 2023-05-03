! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

PROGRAM r4divo_co   ! DO IVO CALC ONLY FOR SMALL BASIS SETS

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module
    use module_file_manager
    use module_2integrals
    use module_realonly, only: check_realonly, realonly
    use read_input_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer                 :: input_unit, nuniq
    logical                 :: test
    character*50            :: filename

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!        debug = .TRUE.
    debug = .FALSE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
#ifdef HAVE_MPI
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_rank(MPI_COMM_WORLD, rank, ierr)
#else
    rank = 0; nprocs = 1
#endif

    tmem = 0.0d+00
    val = 0
    Call DATE_AND_TIME(VALUES=val)
    totalsec = val(8)*(1.0d-03) + val(7) + val(6)*(6.0d+01) + val(5)*(6.0d+01)**2
    initdate = val(3)
    inittime = totalsec
    if (rank == 0) then
        print '(A,I8,A,I8)', 'initialization of mpi, rank :', rank, ' nprocs :', nprocs
        print *, ''
        print *, ' ENTER R4DIVO PROGRAM written by M. Abe 2007/7/20'
        print *, ''
        print '("Current Memory is ",F10.2,"MB")', tmem/1024/1024

        print *, 'Year = ', val(1), 'Mon = ', val(2), 'Date = ', val(3)
        print *, 'Hour = ', val(5), 'Min = ', val(6), 'Sec = ', val(7), '.', val(8)
        print *, 'inittime = ', inittime
    end if

    nhomo = 0  ! Default value of nhomo
    Call timing(val(3), totalsec, date0, tsec)
    call open_formatted_file(unit=input_unit, file='active.inp', status="old", optional_action='read')
    call read_input(input_unit)
    if (rank == 0) then
        print *, 'ninact     =', ninact
        print *, 'nact       =', nact
        print *, 'nsec       =', nsec
        print *, 'nelec      =', nelec
        print *, 'nroot      =', nroot
        print *, 'selectroot =', selectroot
        print *, 'totsym     =', totsym
        print *, 'ncore      =', ncore
        print *, 'nbas       =', nbas
        print *, 'eshift     =', eshift          ! NO USE IN IVO BUT FOR CASCI AND CASPT2 IT IS USED
        print *, 'nhomo      =', nhomo
        print *, 'lscom      =', lscom
        print *, 'noccg      =', noccg
        print *, 'noccu      =', noccu
        print *, 'nvcutg     =', nvcutg
        print *, 'nvcutu     =', nvcutu
        print *, 'diracver   =', dirac_version
    end if

    ! Read MRCONEE file (orbital energies, symmetries and multiplication tables)
    filename = 'MRCONEE'
    call readorb_enesym(filename)
    call read1mo(filename)

    call check_realonly
    ! Create UTChem type MDCINT file from Dirac MDCINT file
    if (rank == 0) print *, "Create_newmdcint"
    call create_newmdcint

    call get_mdcint_filename(0)
    ! Read UTChem type MDCINT files and expands the 2-electron integral in memory
    call readint2_casci(mdcintnew, nuniq)

    realcvec = .TRUE.
    test = .true.
    realc = .FALSE.      !!!      realc =.TRUE.
    realcvec = .FALSE.   !!!      realcvec =.TRUE.
!    This is test for bug fix about realc part
    if (rank == 0) then
        print '("Current Memory is ",F10.2,"MB")', tmem/1024/1024
        print *, ' '
        print *, '*******************************'
        print *, ' '
        print *, 'IREP IS ', repna(totsym)
        print *, ' '
        print *, '*******************************'
        print *, ' '
        print *, realc, 'realc'
        print *, realcvec, 'realcvec'
        print *, 'FOR TEST WE DO (F,F)'
        print *, realc, 'realc'
        print *, realcvec, 'realcvec'
    end if
    iroot = selectroot

    if (rank == 0) then
        Allocate (fock_cmplx(nsec, nsec)); Call memplus(KIND(fock_cmplx), SIZE(fock_cmplx), 2)

        fock_cmplx(:, :) = 0.0d+00

!! NOW MAKE FOCK MATRIX FOR IVO (only virtual spinors)
!! fij = hij + SIGUMA_a(ij|aa)-(ia|aj)}

        Call fockivo

        deallocate (fock_cmplx); Call memminus(KIND(fock_cmplx), SIZE(fock_cmplx), 2)
    end if

    ! Deallocate memory
    deallocate (irpamo); Call memminus(KIND(irpamo), SIZE(irpamo), 1)
    deallocate (indmo_cas_to_dirac); Call memminus(KIND(indmo_cas_to_dirac), SIZE(indmo_cas_to_dirac), 1)
    deallocate (indmo_dirac_to_cas); Call memminus(KIND(indmo_dirac_to_cas), SIZE(indmo_dirac_to_cas), 1)
    deallocate (one_elec_int_i); Call memminus(KIND(one_elec_int_i), SIZE(one_elec_int_i), 1)
    deallocate (one_elec_int_r); Call memminus(KIND(one_elec_int_r), SIZE(one_elec_int_r), 1)
    deallocate (int2r_f1); Call memminus(KIND(int2r_f1), SIZE(int2r_f1), 1)
    deallocate (int2r_f2); Call memminus(KIND(int2r_f2), SIZE(int2r_f2), 1)
    if (.not. realonly%is_realonly()) then
        deallocate (int2i_f1); Call memminus(KIND(int2i_f1), SIZE(int2i_f1), 1)
        deallocate (int2i_f2); Call memminus(KIND(int2i_f2), SIZE(int2i_f2), 1)
    end if
    deallocate (MULTB_S); Call memminus(KIND(MULTB_S), SIZE(MULTB_S), 1)
    deallocate (MULTB_D); Call memminus(KIND(MULTB_D), SIZE(MULTB_D), 1)
    deallocate (MULTB_DS); Call memminus(KIND(MULTB_DS), SIZE(MULTB_DS), 1)
    deallocate (MULTB_DF); Call memminus(KIND(MULTB_DF), SIZE(MULTB_DF), 1)
    deallocate (MULTB_DB); Call memminus(KIND(MULTB_DB), SIZE(MULTB_DB), 1)
    deallocate (MULTB_SB); Call memminus(KIND(MULTB_SB), SIZE(MULTB_SB), 1)

    if (rank == 0) then
        print *, '("Current Memory is ",F10.2,"MB")', tmem/1024/1024
        Call timing(val(3), totalsec, date0, tsec0)
        print *, 'End r4divo part'
    end if
#ifdef HAVE_MPI
    call MPI_FINALIZE(ierr)
#endif

END program r4divo_co
