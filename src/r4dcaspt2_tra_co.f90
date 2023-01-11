! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

PROGRAM r4dcaspt2_tra_co   ! DO CASPT2 CALC WITH MO TRANSFORMATION

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module
    use module_file_manager
    use read_input_module, only: read_input
    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
    real(16)                :: time0, time1
#endif
    integer                 :: ieshift, input_unit = default_unit, new_unit = default_unit
    real*8                  :: e0, e2, e2all, weight0
    complex*16, allocatable :: ci(:)
    real*8, allocatable     :: ecas(:)
    character*50            :: filename
    integer                 :: idetr_array_len ! length of array = idetr(1:2**nact - 1)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!        debug = .TRUE.
    debug = .FALSE.
    thres = 1.0d-15
!        thres = 0.0d+00

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

!   MPI initialization and get the number of MPI processes (nprocs) and own process number.
#ifdef HAVE_MPI
    call MPI_INIT(ierr)
    time0 = MPI_Wtime()
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_rank(MPI_COMM_WORLD, rank, ierr)
#else
    rank = 0; nprocs = 1
#endif
    if (rank == 0) then
        print *, ''
        print *, ' ENTER R4DCASPT2_TRA_TY PROGRAM written by M. Abe 2007.7.23'
        print *, ''
    end if
    tmem = 0.0d+00

    val = 0
    Call DATE_AND_TIME(VALUES=val)
    if (rank == 0) then
        print *, 'Year = ', val(1), 'Mon = ', val(2), 'Date = ', val(3)
        print *, 'Hour = ', val(5), 'Min = ', val(6), 'Sec = ', val(7), '.', val(8)
    end if
    totalsec = val(8)*(1.0d-03) + val(7) + val(6)*(6.0d+01) + val(5)*(6.0d+01)**2
    initdate = val(3)
    inittime = totalsec

    if (rank == 0) then
        print *, inittime
        Call timing(val(3), totalsec, date0, tsec)
    end if

    eshift = 0.0d+00
    ieshift = 0

    call open_formatted_file(unit=input_unit, file='active.inp', status="old", optional_action='read')
    call read_input(input_unit)
    close (input_unit)

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
        if (is_ras1_configured) print *, "RAS1 =", ras1_list
        if (is_ras2_configured) print *, "RAS2 =", ras2_list
        if (is_ras3_configured) print *, "RAS3 =", ras3_list
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    If (mod(nelec, 2) == 0) then
        evenelec = .true.
    else
        evenelec = .false.
    End if

!        print *,' ENTER READ r4dmoint1'
    if (rank == 0) print *, ' ENTER READ MRCONEE'
    filename = 'MRCONEE'

    call readorb_enesym_co(filename)
    call read1mo_co(filename)

    if (rank == 0) then
        print *, ' EXIT READ MRCONEE'
        print *, ' ENTER READ MDCINT'
    end if

    ! Get MDCINTNEWX's filename and subspace filename
    call get_mdcint_filename(0)
    call get_subspace_filename

    Call readint2_ord_co(mdcintnew)

    if (rank == 0) print '("Current Memory is ",F10.2,"MB")', tmem/1024/1024
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    nmo = ninact + nact + nsec
    if (rank == 0) print *, 'nmo        =', nmo

    call open_unformatted_file(unit=new_unit, file="CIMAT", status='old', optional_action="read")

    read (new_unit) ndet
    Allocate (idet(1:ndet)); Call memplus(KIND(idet), SIZE(idet), 1)
    Allocate (ecas(1:ndet)); Call memplus(KIND(ecas), SIZE(ecas), 1)

    read (new_unit) idet(1:ndet)
    read (new_unit) ecas(1:ndet)

    read (new_unit) idetr_array_len
    allocate (idetr(1:idetr_array_len)); call memplus(kind(idet), size(idet), 1)
    read (new_unit) idetr(1:idetr_array_len)
    close (new_unit)

    Allocate (eigen(1:nroot)); Call memplus(KIND(eigen), SIZE(eigen), 1)
    eigen = 0.0d+00
    eigen(1:nroot) = ecas(1:nroot) + ecore

    Deallocate (ecas)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (rank == 0) print *, ' ENTER READ NEWCICOEFF', ndet

    Allocate (ci(1:ndet))
    ci = 0.0d+00

    call open_unformatted_file(unit=new_unit, file="NEWCICOEFF", status='old', optional_action="read")

    read (new_unit) ci(1:ndet)

    close (new_unit)

    Allocate (cir(1:ndet, selectroot:selectroot))
    Allocate (cii(1:ndet, selectroot:selectroot))

    cir(1:ndet, selectroot) = DBLE(ci(1:ndet))
    cii(1:ndet, selectroot) = DIMAG(ci(1:ndet))

    deallocate (ci)

    if (rank == 0) print *, ' EXIT READ NEWCICOEFF'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call open_unformatted_file(unit=new_unit, file="EPS", status='old', optional_action="read")

    read (new_unit) nmo
    Allocate (eps(1:nmo)); Call memplus(KIND(eps), SIZE(eps), 1)
    eps = 0.0d+00
    read (new_unit) eps(1:nmo)

    close (new_unit)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call open_unformatted_file(unit=new_unit, file="TRANSFOCK", status='old', optional_action="read")

    read (new_unit) nmo
    Allocate (f(nmo, nmo)); Call memplus(KIND(f), SIZE(f), 2)
    read (new_unit) f

    close (new_unit)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (rank == 0) then
        print *, ' '
        print *, '*******************************'
        print *, ' '
        print *, 'IREP IS ', repna(totsym)
        print *, ' '
        print *, '*******************************'
        print *, ' '
    end if
    realcvec = .TRUE.

!    This is test for bug fix about realc part

    if (rank == 0) then
        print *, realc, 'realc'
        print *, realcvec, 'realcvec'
    end if
    realc = .FALSE.      !!!      realc =.TRUE.
    realcvec = .FALSE.   !!!      realcvec =.TRUE.

    if (rank == 0) then
        print *, 'FOR TEST WE DO (F,F)'
        print *, realc, 'realc'
        print *, realcvec, 'realcvec'
    end if
!!=============================================!
!                                              !
    iroot = selectroot                       !
!                                              !
!!=============================================!

!      print *,'RECALCULATION OF CASCI ENERGY'
!      Call e0test_v2

    e2 = 0.0d+00

    Call calce0(e0)

    e2all = 0.0d+00

    date1 = initdate
    tsec1 = totalsec

    Call timing(date1, tsec1, date0, tsec0)
    ! date1 = date0
    ! tsec1 = tsec0
    if (rank == 0) print *, 'A1int filename : ', trim(a1int), ' rank', rank

    ! Call intra_3(2, 1, 2, 2, 'A1int')
    Call intra_3(2, 1, 2, 2, a1int)
    if (rank == 0) print *, 'End intra3 A1int'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_3(2, 1, 1, 1, a2int)
    if (rank == 0) print *, 'End intra_3 A2int'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    sumc2local = 0.0d+00
    if (rank == 0) print *, 'Enter solvA'
    Call solvA_ord_ty(e0, e2)
    e2all = e2all + e2
    if (rank == 0) print *, e2all

    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_2(2, 1, 2, 1, bint)
    if (rank == 0) print *, 'End intra_2 Bint'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    sumc2local = 0.0d+00
    Call solvB_ord_ty(e0, e2)
    e2all = e2all + e2
    if (rank == 0) print *, e2all

    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_3(3, 2, 2, 2, c1int)
    if (rank == 0) print *, 'End intra_3 C1int'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_3(3, 2, 1, 1, c2int)
    if (rank == 0) print *, 'End intra_3 C2int'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_1(3, 1, 1, 2, c3int)
    if (rank == 0) print *, 'End intra_1 C3int'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    sumc2local = 0.0d+00
    Call solvC_ord_ty(e0, e2)
    e2all = e2all + e2

    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_3(3, 1, 2, 2, d1int)
    if (rank == 0) print *, 'End intra_1 D1int'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_1(3, 2, 2, 1, d2int)
    if (rank == 0) print *, 'End intra_1 D2int'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_3(3, 1, 1, 1, d3int)
    if (rank == 0) print *, 'End intra_1 D3int'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    sumc2local = 0.0d+00
    Call solvD_ord_ty(e0, e2)
    e2all = e2all + e2
    if (rank == 0) print *, e2all

    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_1(3, 1, 2, 1, eint)
    if (rank == 0) print *, 'End intra_1 Eint'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    sumc2local = 0.0d+00
    Call solvE_ord_ty(e0, e2)
    e2all = e2all + e2
    if (rank == 0) print *, e2all

    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_2(3, 2, 3, 2, fint)
    if (rank == 0) print *, 'End intra_1 Fint'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    sumc2local = 0.0d+00
    Call solvF_ord_ty(e0, e2)
    e2all = e2all + e2
    if (rank == 0) print *, e2all

    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_1(3, 1, 3, 2, gint)
    if (rank == 0) print *, 'End intra_1 Gint'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    sumc2local = 0.0d+00
    Call solvG_ord_ty(e0, e2)
    e2all = e2all + e2
    if (rank == 0) print *, e2all
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    sumc2local = 0.0d+00
    if (rank == 0) print *, 'enter solveH_ord_ty'
    if (rank == 0) print *, 'Enter intra_2 Hint'
    Call intra_2(3, 1, 3, 1, hint)
    if (rank == 0) print *, 'End intra_2 Hint'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call solvH_ord_ty(e0, e2)
    e2all = e2all + e2
    if (rank == 0) print *, e2all

    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    if (rank == 0) print '("c^2 ",F30.15)', sumc2
    weight0 = 1.0d+00/(1.0d+00 + sumc2)

    if (rank == 0) then
        print '("weight of 0th wave function is",F30.15)', weight0

        print '("Total second order energy is ",F30.15," a.u.")', e2all - eshift*sumc2
        print '("Total energy is ",F30.15," a.u.")', e2all + eigen(iroot) - eshift*sumc2
    end if
    if (allocated(cir)) deallocate (cir); Call memminus(KIND(cir), SIZE(cir), 1)
    if (allocated(cii)) deallocate (cii); Call memminus(KIND(cii), SIZE(cii), 1)
    if (allocated(eigen)) deallocate (eigen); Call memminus(KIND(eigen), SIZE(eigen), 1)
    if (allocated(eps)) deallocate (eps); Call memminus(KIND(eps), SIZE(eps), 1)
    if (allocated(idet)) deallocate (idet); Call memminus(KIND(idet), SIZE(idet), 1)
    if (allocated(idetr)) deallocate (idetr); Call memminus(KIND(MULTB_S), SIZE(MULTB_S), 1)
    if (allocated(MULTB_S)) deallocate (MULTB_S); Call memminus(KIND(MULTB_S), SIZE(MULTB_S), 1)
    if (allocated(MULTB_D)) deallocate (MULTB_D); Call memminus(KIND(MULTB_D), SIZE(MULTB_D), 1)
    if (allocated(MULTB_DS)) deallocate (MULTB_DS); Call memminus(KIND(MULTB_DS), SIZE(MULTB_DS), 1)
    if (allocated(MULTB_DF)) deallocate (MULTB_DF); Call memminus(KIND(MULTB_DF), SIZE(MULTB_DF), 1)
    if (allocated(MULTB_DB)) deallocate (MULTB_DB); Call memminus(KIND(MULTB_DB), SIZE(MULTB_DB), 1)
    if (allocated(MULTB_SB)) deallocate (MULTB_SB); Call memminus(KIND(MULTB_SB), SIZE(MULTB_SB), 1)

    Call timing(val(3), totalsec, date0, tsec0)
    if (rank == 0) print *, 'End r4dcaspt2_tra_ty'
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    time1 = MPI_Wtime()
    if (rank == 0) then
        write (*, "(a,e16.6)") "MPI_Wtime :", time1 - time0
    end if
    call MPI_FINALIZE(ierr)
#endif

END program r4dcaspt2_tra_co
