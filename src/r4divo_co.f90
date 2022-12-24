! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

PROGRAM r4divo_co   ! DO IVO CALC ONLY FOR SMALL BASIS SETS

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module
    use module_file_manager
    use read_input_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer                 :: input_unit
    logical                 :: test
    character*50            :: filename

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!        debug = .TRUE.
    debug = .FALSE.
    thres = 1.0d-15
!        thres = 0.0d+00

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
        print *, ' ENTER R4DIVO PROGRAM written by M. Abe test17_ty version 2007/7/20'
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
        print *, 'diracver   =', dirac_version
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    filename = 'MRCONEE'

    call readorb_enesym_co(filename)

    call read1mo_co(filename)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Iwamuro create new ikr for dirac
    if (rank == 0) print *, "Create_newmdcint"
    call create_newmdcint

    call get_mdcint_filename(0)
    Call readint2_ivo_co(mdcintnew)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!=============================================!
!                                              !
    iroot = selectroot
!                                              !
!!=============================================!

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!            BUILDING  FOCK MATRIX               !
!  fij = hij + SIGUMA[<0|Ekl|0>{(ij|kl)-(il|kj)} !
!                 kl                             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

!! TEST TO CALCULATE FOCK MATRIX OF HF STATE fpq = hpq + SIGUMA_r[(pq|rr)-(pr|qr)]
!! THIS MUST BE DIAGONAL MATRIX AND DIAGONAL ELEMENTS CORESPONDS TO SPINOR ENERGIES.
    if (rank == 0) then
        Allocate (f(nsec, nsec)); Call memplus(KIND(f), SIZE(f), 2)

        f(:, :) = 0.0d+00

!! NOW MAKE FOCK MATRIX FOR IVO (only virtual spinors

!! fij = hij + SIGUMA_a(ij|aa)-(ia|aj)}

        Call fockivo_co

        deallocate (f); Call memminus(KIND(f), SIZE(f), 2)
    end if
    if (allocated(MULTB_S)) deallocate (MULTB_S); Call memminus(KIND(MULTB_S), SIZE(MULTB_S), 1)
    if (allocated(MULTB_D)) deallocate (MULTB_D); Call memminus(KIND(MULTB_D), SIZE(MULTB_D), 1)
    if (allocated(MULTB_DS)) deallocate (MULTB_DS); Call memminus(KIND(MULTB_DS), SIZE(MULTB_DS), 1)
    if (allocated(MULTB_DF)) deallocate (MULTB_DF); Call memminus(KIND(MULTB_DF), SIZE(MULTB_DF), 1)
    if (allocated(MULTB_DB)) deallocate (MULTB_DB); Call memminus(KIND(MULTB_DB), SIZE(MULTB_DB), 1)
    if (allocated(MULTB_SB)) deallocate (MULTB_SB); Call memminus(KIND(MULTB_SB), SIZE(MULTB_SB), 1)
    if (allocated(irpmo)) deallocate (irpmo); Call memminus(KIND(irpmo), SIZE(irpmo), 1)
    if (allocated(irpamo)) deallocate (irpamo); Call memminus(KIND(irpamo), SIZE(irpamo), 1)
    if (allocated(indmo)) deallocate (indmo); Call memminus(KIND(indmo), SIZE(indmo), 1)
    if (allocated(indmor)) deallocate (indmor); Call memminus(KIND(indmor), SIZE(indmor), 1)
    if (allocated(onei)) deallocate (onei); Call memminus(KIND(onei), SIZE(onei), 1)
    if (allocated(oner)) deallocate (oner); Call memminus(KIND(oner), SIZE(oner), 1)
    if (allocated(int2r_f1)) deallocate (int2r_f1); Call memminus(KIND(int2r_f1), SIZE(int2r_f1), 1)
    if (allocated(int2i_f1)) deallocate (int2i_f1); Call memminus(KIND(int2i_f1), SIZE(int2i_f1), 1)
    if (allocated(int2r_f2)) deallocate (int2r_f2); Call memminus(KIND(int2r_f2), SIZE(int2r_f2), 1)
    if (allocated(int2i_f2)) deallocate (int2i_f2); Call memminus(KIND(int2i_f2), SIZE(int2i_f2), 1)

    if (rank == 0) then
        print '("Current Memory is ",F10.2,"MB")', tmem/1024/1024

        Call timing(val(3), totalsec, date0, tsec0)
        print *, 'End r4divo_co part'
    end if
#ifdef HAVE_MPI
    call MPI_FINALIZE(ierr)
#endif

END program r4divo_co
