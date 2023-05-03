! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

PROGRAM r4dcaspt2_tra   ! DO CASPT2 CALC WITH MO TRANSFORMATION

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module
    use module_dict, only: add
    use module_file_manager
    use module_realonly, only: check_realonly, realonly
    use read_input_module, only: read_input
    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
    real(16)                :: time0, time1
#endif
    integer                 :: ieshift, unit_input, unit_new
    real*8                  :: e0, e2, e2all, weight0
    complex*16, allocatable :: ci(:)
    real*8, allocatable     :: ecas(:)
    character*50            :: filename
    integer                 :: dict_size ! The number of CAS configurations
    integer                 :: idx, dict_key, dict_val

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    debug = .FALSE.

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
        print '(2(A,1X,I0))', 'initialization of mpi, rank :', rank, ' nprocs :', nprocs
        print *, ''
        print *, ' ENTER R4DCASPT2_TRA PROGRAM written by M. Abe 2007.7.23'
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

    If (mod(nelec, 2) == 0) then
        evenelec = .true.
    else
        evenelec = .false.
    End if

    ! Read MRCONEE file (orbital energies, symmetries and multiplication tables)
    if (rank == 0) print *, ' ENTER READ MRCONEE'
    filename = 'MRCONEE'
    call readorb_enesym(filename)
    call read1mo(filename)
    call check_realonly()

    if (rank == 0) then
        print *, ' EXIT READ MRCONEE'
        print *, ' ENTER READ MDCINT'
    end if

    ! Get MDCINTNEWx's filename and subspace filename
    call get_mdcint_filename(0)
    call get_subspace_filename

    ! READ MDCINTNEWx's file and devide into each subspace files.
    Call readint2_ord_co(mdcintnew)

    if (rank == 0) print '("Current Memory is ",F10.2,"MB")', tmem/1024/1024
    nmo = ninact + nact + nsec
    if (rank == 0) print *, 'nmo        =', nmo

    ! Read CAS configuration convertion list
    call open_unformatted_file(unit=unit_new, file="CIMAT", status='old', optional_action="read")
    read (unit_new) ndet
    Allocate (cas_idx(1:ndet)); Call memplus(KIND(cas_idx), SIZE(cas_idx), 1)
    Allocate (ecas(1:ndet)); Call memplus(KIND(ecas), SIZE(ecas), 1)
    read (unit_new) cas_idx(1:ndet)
    read (unit_new) ecas(1:ndet)
    read (unit_new) dict_size ! The number of CAS configurations
    do idx = 1, dict_size
        read (unit_new) dict_key, dict_val
        call add(dict_cas_idx_reverse, dict_key, dict_val)
    end do
    close (unit_new)

    ! Read CASCI energy
    Allocate (eigen(1:nroot)); Call memplus(KIND(eigen), SIZE(eigen), 1)
    eigen = 0.0d+00
    eigen(1:nroot) = ecas(1:nroot) + ecore
    Deallocate (ecas)

    ! Read CI coefficients
    if (rank == 0) print *, ' ENTER READ NEWCICOEFF', ndet
    Allocate (ci(1:ndet))
    ci = 0.0d+00
    call open_unformatted_file(unit=unit_new, file="NEWCICOEFF", status='old', optional_action="read")
    read (unit_new) ci(1:ndet)
    close (unit_new)
    Allocate (cir(1:ndet, selectroot:selectroot))
    Allocate (cii(1:ndet, selectroot:selectroot))
    cir(1:ndet, selectroot) = DBLE(ci(1:ndet))
    cii(1:ndet, selectroot) = DIMAG(ci(1:ndet))
    deallocate (ci)
    if (rank == 0) print *, ' EXIT READ NEWCICOEFF'

    ! Read epsilons
    call open_unformatted_file(unit=unit_new, file="EPS", status='old', optional_action="read")
    read (unit_new) nmo
    Allocate (eps(1:nmo)); Call memplus(KIND(eps), SIZE(eps), 1)
    eps = 0.0d+00
    read (unit_new) eps(1:nmo)
    close (unit_new)
    ! Read the transformed Fock matrix
    call open_unformatted_file(unit=unit_new, file="TRANSFOCK", status='old', optional_action="read")
    read (unit_new) nmo
    if (realonly%is_realonly()) then
        allocate (fock_real(nmo, nmo)); Call memplus(KIND(fock_real), SIZE(fock_real), 1)
        read (unit_new) fock_real
        allocate (fock_cmplx(nmo, nmo)); Call memplus(KIND(fock_cmplx), SIZE(fock_cmplx), 2)
        fock_cmplx = 0.0d+00
        fock_cmplx = fock_real
    else
        Allocate (fock_cmplx(nmo, nmo)); Call memplus(KIND(fock_cmplx), SIZE(fock_cmplx), 2)
        read (unit_new) fock_cmplx
    end if

    close (unit_new)

    ! Print the irreducible representation that calculates energy
    if (rank == 0) then
        print *, ' '
        print *, '*******************************'
        print *, ' '
        print *, 'IREP IS ', repna(totsym)
        print *, ' '
        print *, '*******************************'
        print *, ' '
    end if
    iroot = selectroot
    realc = .FALSE.      !!!      realc =.TRUE.
    realcvec = .FALSE.   !!!      realcvec =.TRUE.

    ! Calculate the 0th order energy
    e2 = 0.0d+00
    Call calce0(e0)

    ! Initialize the date, time and the 2nd order energy
    e2all = 0.0d+00
    date1 = initdate
    tsec1 = totalsec
    Call timing(date1, tsec1, date0, tsec0)

    ! Transform A subspace 2-electron integrals (active, inactive | active, active)
    if (rank == 0) print *, 'A1int filename : ', trim(a1int), ' rank', rank
    Call intra_3(2, 1, 2, 2, a1int)
    if (rank == 0) print *, 'End intra3 A1int'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Transform A subspace 2-electron integrals (active, inactive | inactive, inactive)
    Call intra_3(2, 1, 1, 1, a2int)
    if (rank == 0) print *, 'End intra_3 A2int'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Calculate the A subspace 2nd order energy
    sumc2local = 0.0d+00
    if (rank == 0) print *, 'Enter solvA'
    Call solve_A_subspace(e0, e2)
    e2all = e2all + e2
    if (rank == 0) print *, e2all
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Transform B subspace 2-electron integrals (active, inactive | active, inactive)
    Call intra_2(2, 1, 2, 1, bint)
    if (rank == 0) print *, 'End intra_2 Bint'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Calculate the B subspace 2nd order energy
    sumc2local = 0.0d+00
    Call solve_B_subspace(e0, e2)
    e2all = e2all + e2
    if (rank == 0) print *, e2all
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Transform C subspace 2-electron integrals (secondary, active | active, active)
    Call intra_3(3, 2, 2, 2, c1int)
    if (rank == 0) print *, 'End intra_3 C1int'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Transform C subspace 2-electron integrals (secondary, active | inactive, inactive)
    Call intra_3(3, 2, 1, 1, c2int)
    if (rank == 0) print *, 'End intra_3 C2int'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Transform C subspace 2-electron integrals (secondary, inactive | inactive, active)
    Call intra_1(3, 1, 1, 2, c3int)
    if (rank == 0) print *, 'End intra_1 C3int'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Calculate the B subspace 2nd order energy
    sumc2local = 0.0d+00
    Call solve_C_subspace(e0, e2)
    e2all = e2all + e2
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Transform D subspace 2-electron integrals (secondary, inactive | active, active)
    Call intra_3(3, 1, 2, 2, d1int)
    if (rank == 0) print *, 'End intra_1 D1int'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Transform D subspace 2-electron integrals (secondary, active | active, inactive)
    Call intra_1(3, 2, 2, 1, d2int)
    if (rank == 0) print *, 'End intra_1 D2int'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Transform D subspace 2-electron integrals (secondary, inactive | inactive, inactive)
    Call intra_3(3, 1, 1, 1, d3int)
    if (rank == 0) print *, 'End intra_1 D3int'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Calculate the D subspace 2nd order energy
    sumc2local = 0.0d+00
    Call solve_D_subspace(e0, e2)
    e2all = e2all + e2
    if (rank == 0) print *, e2all
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Transform E subspace 2-electron integrals (secondary, active | active, inactive)
    Call intra_1(3, 1, 2, 1, eint)
    if (rank == 0) print *, 'End intra_1 Eint'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Calculate the E subspace 2nd order energy
    sumc2local = 0.0d+00
    Call solve_E_subspace(e0, e2)
    e2all = e2all + e2
    if (rank == 0) print *, e2all
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Transform F subspace 2-electron integrals (secondary, active | secondary, active)
    Call intra_2(3, 2, 3, 2, fint)
    if (rank == 0) print *, 'End intra_1 Fint'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Calculate the F subspace 2nd order energy
    sumc2local = 0.0d+00
    Call solve_F_subspace(e0, e2)
    e2all = e2all + e2
    if (rank == 0) print *, e2all
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Transform G subspace 2-electron integrals (secondary, inactive | secondary, active)
    Call intra_1(3, 1, 3, 2, gint)
    if (rank == 0) print *, 'End intra_1 Gint'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Calculate the G subspace 2nd order energy
    sumc2local = 0.0d+00
    Call solve_G_subspace(e0, e2)
    e2all = e2all + e2
    if (rank == 0) print *, e2all
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Transform H subspace 2-electron integrals (secondary, inactive | secondary, inactive)
    if (rank == 0) print *, 'Enter intra_2 Hint'
    Call intra_2(3, 1, 3, 1, hint)
    if (rank == 0) print *, 'End intra_2 Hint'
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Calculate the H subspace 2nd order energy
    sumc2local = 0.0d+00
    Call solve_H_subspace(e0, e2)
    e2all = e2all + e2
    if (rank == 0) print *, e2all
    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    ! Print out the total 2nd order energy
    if (rank == 0) print '("c^2 ",F30.15)', sumc2
    ! Calculate and print the weight of the 0th wave function
    weight0 = 1.0d+00/(1.0d+00 + sumc2)
    if (rank == 0) then
        print '("weight of 0th wave function is",F30.15)', weight0

        print '("Total second order energy is ",F30.15," a.u.")', e2all - eshift*sumc2
        print '("Total energy is ",F30.15," a.u.")', e2all + eigen(iroot) - eshift*sumc2
    end if

    ! Deallocate the memory
    if (allocated(cir)) Call memminus(KIND(cir), SIZE(cir), 1); deallocate (cir)
    if (allocated(cii)) Call memminus(KIND(cii), SIZE(cii), 1); deallocate (cii)
    if (allocated(eigen)) Call memminus(KIND(eigen), SIZE(eigen), 1); deallocate (eigen)
    if (allocated(eps)) Call memminus(KIND(eps), SIZE(eps), 1); deallocate (eps)
    if (allocated(cas_idx)) Call memminus(KIND(cas_idx), SIZE(cas_idx), 1); deallocate (cas_idx)
    if (allocated(MULTB_S)) Call memminus(KIND(MULTB_S), SIZE(MULTB_S), 1); deallocate (MULTB_S)
    if (allocated(MULTB_D)) Call memminus(KIND(MULTB_D), SIZE(MULTB_D), 1); deallocate (MULTB_D)
    if (allocated(MULTB_DS)) Call memminus(KIND(MULTB_DS), SIZE(MULTB_DS), 1); deallocate (MULTB_DS)
    if (allocated(MULTB_DF)) Call memminus(KIND(MULTB_DF), SIZE(MULTB_DF), 1); deallocate (MULTB_DF)
    if (allocated(MULTB_DB)) Call memminus(KIND(MULTB_DB), SIZE(MULTB_DB), 1); deallocate (MULTB_DB)
    if (allocated(MULTB_SB)) Call memminus(KIND(MULTB_SB), SIZE(MULTB_SB), 1); deallocate (MULTB_SB)

    ! Print out the total time
    Call timing(val(3), totalsec, date0, tsec0)
    if (rank == 0) print *, 'End r4dcaspt2_tra'
#ifdef HAVE_MPI
    ! MPI finalization
    time1 = MPI_Wtime()
    if (rank == 0) then
        ! Print out the total time (MPI)
        write (*, "(a,e16.6)") "MPI_Wtime :", time1 - time0
    end if
    call MPI_FINALIZE(ierr)
#endif

END program r4dcaspt2_tra
