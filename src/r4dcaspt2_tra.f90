! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

PROGRAM r4dcaspt2_tra   ! DO CASPT2 CALC WITH MO TRANSFORMATION

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_dict, only: add
    use module_error, only: stop_with_errorcode
    use module_file_manager
    use module_global_variables
    use module_intra, only: intra_1, intra_2, intra_3
    use module_realonly, only: check_realonly, realonly
    use module_time
    use read_input_module, only: read_input
    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
    real(16)                :: time0, time1
#endif
    integer                 :: unit_input, unit_new
    real(8)                 :: e0, e2, e2all, weight0
    complex*16, allocatable         :: ci(:)
    real(8), allocatable            :: ecas(:)
    character(:), allocatable       :: filename
    integer                 :: dict_cas_idx_size, dict_cas_idx_reverse_size ! The number of CAS configurations
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

    call write_allocated_memory_size
    call get_current_time(init_time); call print_time(init_time); start_time = init_time

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
        print *, 'mdcint_scheme =', mdcint_scheme
        if (ras1_size /= 0) print *, "RAS1 =", ras1_list
        if (ras2_size /= 0) print *, "RAS2 =", ras2_list
        if (ras3_size /= 0) print *, "RAS3 =", ras3_list
    end if

    if (ninact == 0 .and. nsec == 0) then
        if (rank == 0) print *, "The CASPT2 energy cannot be defined when ninact = 0 and nsec = 0."
        stop
    end if

    ! Read MRCONEE file (orbital energies, symmetries and multiplication tables)
    if (rank == 0) print *, ' ENTER READ MRCONEE'
    filename = 'MRCONEE'
    call check_dirac_integer_size(filename)
    call read_mrconee(filename)
    call check_realonly()

    if (rank == 0) then
        print *, ' EXIT READ MRCONEE'
        print *, ' ENTER READ MDCINT'
    end if

    ! Get MDCINTNEWx's filename and subspace filename
    call get_mdcint_filename(0)
    call get_subspace_filename

    ! READ MDCINTNEWx's file and devide into each subspace files.
    Call divide_2_elec_integral_into_subspaces(mdcintnew)

    call write_allocated_memory_size
    nmo = ninact + nact + nsec
    if (rank == 0) print *, 'nmo        =', nmo

    ! Read CAS configuration convertion list
    call open_unformatted_file(unit=unit_new, file="CIMAT", status='old', optional_action="read")
    read (unit_new) ndet
    Allocate (ecas(1:ndet)); Call memplus(KIND(ecas), SIZE(ecas), 1)
    read (unit_new) ecas(1:ndet)
    read (unit_new) dict_cas_idx_size ! The number of CAS configurations
    do idx = 1, dict_cas_idx_size
        read (unit_new) dict_key, dict_val
        call add(dict_cas_idx, dict_key, dict_val)
    end do
    read (unit_new) dict_cas_idx_reverse_size ! The number of CAS configurations
    do idx = 1, dict_cas_idx_reverse_size
        read (unit_new) dict_key, dict_val
        call add(dict_cas_idx_reverse, dict_key, dict_val)
    end do
    close (unit_new)
    ! Check if dict_cas_idx_size is equal to ndet
    if (dict_cas_idx_size /= ndet .or. dict_cas_idx_reverse_size /= ndet) then
        if (rank == 0) print *, 'ERROR: dict_cas_idx_size /= ndet .or. dict_cas_idx_reverse_size /= ndet. ndet =', ndet, &
            ",dict_cas_idx_size =", dict_cas_idx_size, ",dict_cas_idx_reverse_size =", dict_cas_idx_reverse_size
        call stop_with_errorcode(1)
    end if

    ! Read CASCI energy
    Allocate (eigen(1:nroot)); Call memplus(KIND(eigen), SIZE(eigen), 1)
    eigen = 0.0d+00
    eigen(1:nroot) = ecas(1:nroot) + ecore
    Deallocate (ecas)

    ! Read CI coefficients
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

    ! Calculate eigenvalues of a 0th-order Hamiltonian applied to a 0th-order wave function
    Call calce0(e0)

    ! Initialize the date, time and the 2nd order energy
    e2 = 0.0d+00
    e2all = 0.0d+00
    call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

    if (ninact == 0) then
        if (rank == 0) print *, "Skip the calculation of A subspace 2nd order energy &
&        because the 2nd order energy of A subspace cannot be defined when ninact = 0."
    else
        ! Transform A subspace 2-electron integrals (active, inactive | active, active)
        Call intra_3(2, 1, 2, 2, a1int)
        if (rank == 0) print *, 'End intra3 A1int'
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Transform A subspace 2-electron integrals (active, inactive | inactive, inactive)
        Call intra_3(2, 1, 1, 1, a2int)
        if (rank == 0) print *, 'End intra_3 A2int'
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Calculate the A subspace 2nd order energy
        sumc2local = 0.0d+00
        if (rank == 0) print *, "Start calcultion of A subspace 2nd order energy"
        Call solve_A_subspace(e0, e2)
        e2all = e2all + e2
        if (rank == 0) print *, "End calcultion of A subspace 2nd order energy"
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
    end if

    if (ninact == 0) then
        if (rank == 0) print *, "Skip the calculation of B subspace 2nd order energy &
&        because the 2nd order energy of B subspace cannot be defined when ninact = 0."
    else
        ! Transform B subspace 2-electron integrals (active, inactive | active, inactive)
        Call intra_2(2, 1, 2, 1, bint)
        if (rank == 0) print *, "End intra_2 Bint"
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Calculate the B subspace 2nd order energy
        sumc2local = 0.0d+00
        if (rank == 0) print *, "Start calcuation of B subspace 2nd order energy"
        Call solve_B_subspace(e0, e2)
        e2all = e2all + e2
        if (rank == 0) print *, "End calcuation of B subspace 2nd order energy"
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
    end if

    if (nsec == 0) then
        if (rank == 0) print *, "Skip the calculation of C subspace 2nd order energy &
&        because the 2nd order energy of C subspace cannot be defined when nsec = 0."
    else
        ! Transform C subspace 2-electron integrals (secondary, active | active, active)
        Call intra_3(3, 2, 2, 2, c1int)
        if (rank == 0) print *, 'End intra_3 C1int'
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Transform C subspace 2-electron integrals (secondary, active | inactive, inactive)
        Call intra_3(3, 2, 1, 1, c2int)
        if (rank == 0) print *, 'End intra_3 C2int'
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Transform C subspace 2-electron integrals (secondary, inactive | inactive, active)
        Call intra_1(3, 1, 1, 2, c3int)
        if (rank == 0) print *, 'End intra_1 C3int'
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Calculate the B subspace 2nd order energy
        sumc2local = 0.0d+00
        if (rank == 0) print *, "Start calcuation of C subspace 2nd order energy"
        Call solve_C_subspace(e0, e2)
        e2all = e2all + e2
        if (rank == 0) print *, "End calcuation of C subspace 2nd order energy"
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
    end if

    if (ninact == 0 .or. nsec == 0) then
        if (rank == 0) print *, "Skip the calculation of D subspace 2nd order energy &
&        because the 2nd order energy of D subspace cannot be defined when ninact = 0 or nsec = 0."
    else
        ! Transform D subspace 2-electron integrals (secondary, inactive | active, active)
        Call intra_3(3, 1, 2, 2, d1int)
        if (rank == 0) print *, 'End intra_1 D1int'
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Transform D subspace 2-electron integrals (secondary, active | active, inactive)
        Call intra_1(3, 2, 2, 1, d2int)
        if (rank == 0) print *, 'End intra_1 D2int'
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Transform D subspace 2-electron integrals (secondary, inactive | inactive, inactive)
        Call intra_3(3, 1, 1, 1, d3int)
        if (rank == 0) print *, 'End intra_1 D3int'
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Calculate the D subspace 2nd order energy
        sumc2local = 0.0d+00
        if (rank == 0) print *, "Start calcuation of D subspace 2nd order energy"
        Call solve_D_subspace(e0, e2)
        e2all = e2all + e2
        if (rank == 0) print *, "End calcuation of D subspace 2nd order energy"
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
    end if

    if (ninact == 0 .or. nsec == 0) then
        if (rank == 0) print *, "Skip the calculation of E subspace 2nd order energy &
&        because the 2nd order energy of E subspace cannot be defined when ninact = 0 or nsec = 0."
    else
        ! Transform E subspace 2-electron integrals (secondary, active | active, inactive)
        Call intra_1(3, 1, 2, 1, eint)
        if (rank == 0) print *, 'End intra_1 Eint'
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Calculate the E subspace 2nd order energy
        sumc2local = 0.0d+00
        if (rank == 0) print *, "Start calcuation of E subspace 2nd order energy"
        Call solve_E_subspace(e0, e2)
        e2all = e2all + e2
        if (rank == 0) print *, "End calcuation of E subspace 2nd order energy"
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
    end if

    if (nsec == 0) then
        if (rank == 0) print *, "Skip the calculation of F subspace 2nd order energy &
&        because the 2nd order energy of F subspace cannot be defined when nsec = 0."
    else
        ! Transform F subspace 2-electron integrals (secondary, active | secondary, active)
        Call intra_2(3, 2, 3, 2, fint)
        if (rank == 0) print *, 'End intra_1 Fint'
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Calculate the F subspace 2nd order energy
        sumc2local = 0.0d+00
        if (rank == 0) print *, "Start calcuation of F subspace 2nd order energy"
        Call solve_F_subspace(e0, e2)
        e2all = e2all + e2
        if (rank == 0) print *, "End calcuation of F subspace 2nd order energy"
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
    end if

    if (ninact == 0 .or. nsec == 0) then
        if (rank == 0) print *, "Skip the calculation of G subspace 2nd order energy &
&        because the 2nd order energy of G subspace cannot be defined when ninact = 0 or nsec = 0."
    else
        ! Transform G subspace 2-electron integrals (secondary, inactive | secondary, active)
        Call intra_1(3, 1, 3, 2, gint)
        if (rank == 0) print *, 'End intra_1 Gint'
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Calculate the G subspace 2nd order energy
        sumc2local = 0.0d+00
        if (rank == 0) print *, "Start calcuation of G subspace 2nd order energy"
        Call solve_G_subspace(e0, e2)
        e2all = e2all + e2
        if (rank == 0) print *, "End calcuation of G subspace 2nd order energy"
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
    end if

    if (ninact == 0 .or. nsec == 0) then
        if (rank == 0) print *, "Skip the calculation of H subspace 2nd order energy &
&        because the 2nd order energy of H subspace cannot be defined when ninact = 0 or nsec = 0."
    else
        ! Transform H subspace 2-electron integrals (secondary, inactive | secondary, inactive)
        if (rank == 0) print *, 'Enter intra_2 Hint'
        Call intra_2(3, 1, 3, 1, hint)
        if (rank == 0) print *, 'End intra_2 Hint'
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Calculate the H subspace 2nd order energy
        sumc2local = 0.0d+00
        if (rank == 0) print *, "Start calcuation of H subspace 2nd order energy"
        Call solve_H_subspace(e0, e2)
        e2all = e2all + e2
        if (rank == 0) print *, "End calcuation of H subspace 2nd order energy"
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
    end if
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
    if (allocated(MULTB_S)) Call memminus(KIND(MULTB_S), SIZE(MULTB_S), 1); deallocate (MULTB_S)
    if (allocated(MULTB_D)) Call memminus(KIND(MULTB_D), SIZE(MULTB_D), 1); deallocate (MULTB_D)
    if (allocated(MULTB_DS)) Call memminus(KIND(MULTB_DS), SIZE(MULTB_DS), 1); deallocate (MULTB_DS)

    ! Print out the total time
    call get_current_time_and_print_diff(init_time, end_time)
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
