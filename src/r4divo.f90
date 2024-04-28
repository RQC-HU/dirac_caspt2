! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

PROGRAM r4divo_co   ! DO IVO CALC ONLY FOR SMALL BASIS SETS

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables
    use module_file_manager
    use module_2integrals
    use module_realonly, only: check_realonly, realonly
    use module_time
    use read_input_module
    use module_ivo_consistency_check

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer                     :: input_unit, nuniq
    character(:), allocatable   :: filename

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

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

    if (rank == 0) then
        call print_head
        print *, ''
        print *, 'START RELATIVISTIC IVO PROGRAM'
    end if

    tmem = 0.0d+00
    call write_allocated_memory_size
    call get_current_time(init_time); call print_time(init_time); start_time = init_time

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
        print *, 'eshift     =', eshift          ! NO USE IN IVO BUT FOR CASCI AND CASPT2 IT IS USED
        print *, 'nhomo      =', nhomo
        if (inversion) then
            print *, "noccg      =", occ_mo_num(1)
            print *, "noccu      =", occ_mo_num(2)
            print *, "nvcutg     =", vcut_mo_num(1)
            print *, "nvcutu     =", vcut_mo_num(2)
        else
            print *, "nocc      =", occ_mo_num(1)
            print *, "nvcut     =", vcut_mo_num(1)
        end if
        print *, 'diracver   =', dirac_version
        print *, 'scheme     =', mdcint_scheme
        print *, 'debugprint    =', debug
    end if

    ! Read MRCONEE file (orbital energies, symmetries and multiplication tables)
    filename = 'MRCONEE'
    call read_mrconee(filename)

    ! Check consistency of IVO input and DFPCMO file.
    call ivo_consistency_check

    call check_realonly
    ! Create UTChem type MDCINT file from Dirac MDCINT file
    if (rank == 0) print *, "Create_newmdcint"
    call create_newmdcint

    call get_mdcint_filename(0)
    ! Read UTChem type MDCINT files and expands the 2-electron integral in memory
    call readint2_casci(mdcintnew, nuniq)

    if (rank == 0) then
        call write_allocated_memory_size
        print *, ' '
        print *, '*******************************'
        print *, ' '
        print *, 'IREP IS ', repna(totsym)
        print *, ' '
        print *, '*******************************'
    end if
    iroot = selectroot

    if (rank == 0) then
        Allocate (fock_cmplx(nsec, nsec)); Call memplus(KIND(fock_cmplx), SIZE(fock_cmplx), 2)

        fock_cmplx(:, :) = 0.0d+00

!! NOW MAKE FOCK MATRIX FOR IVO (only virtual spinors)
!! fij = hij + SIGUMA_a(ij|aa)-(ia|aj)}

        Call fockivo

        Call memminus(KIND(fock_cmplx), SIZE(fock_cmplx), 2); deallocate (fock_cmplx)
    end if

    ! Deallocate memory
    Call memminus(KIND(irpamo), SIZE(irpamo), 1); deallocate (irpamo)
    Call memminus(KIND(indmo_cas_to_dirac), SIZE(indmo_cas_to_dirac), 1); deallocate (indmo_cas_to_dirac)
    Call memminus(KIND(indmo_dirac_to_cas), SIZE(indmo_dirac_to_cas), 1); deallocate (indmo_dirac_to_cas)
    Call memminus(KIND(one_elec_int_i), SIZE(one_elec_int_i), 1); deallocate (one_elec_int_i)
    Call memminus(KIND(one_elec_int_r), SIZE(one_elec_int_r), 1); deallocate (one_elec_int_r)
    Call memminus(KIND(int2r_f1), SIZE(int2r_f1), 1); deallocate (int2r_f1)
    Call memminus(KIND(int2r_f2), SIZE(int2r_f2), 1); deallocate (int2r_f2)
    if (.not. realonly%is_realonly()) then
        Call memminus(KIND(int2i_f1), SIZE(int2i_f1), 1); deallocate (int2i_f1)
        Call memminus(KIND(int2i_f2), SIZE(int2i_f2), 1); deallocate (int2i_f2)
    end if
    Call memminus(KIND(MULTB_S), SIZE(MULTB_S), 1); deallocate (MULTB_S)
    Call memminus(KIND(MULTB_D), SIZE(MULTB_D), 1); deallocate (MULTB_D)
    Call memminus(KIND(MULTB_DS), SIZE(MULTB_DS), 1); deallocate (MULTB_DS)

    call write_allocated_memory_size
    call get_current_time_and_print_diff(init_time, end_time) ! print the total time
    if (rank == 0) print *, 'END OF RELATIVISTIC IVO PROGRAM'
#ifdef HAVE_MPI
    call MPI_FINALIZE(ierr)
#endif

END program r4divo_co
