! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

subroutine r4dcasci   ! DO CASCI CALC IN THIS PROGRAM!

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use dcaspt2_restart_file, only: read_and_validate_restart_file
    use module_global_variables
    use module_file_manager
    use module_2integrals
    use module_realonly
    use module_time
    use read_input_module

    Implicit NONE
    integer                    :: i0, nuniq, unit_eps, unit_input
    character(:), allocatable  :: filename
    character(*), parameter    :: int_input_form = '(1x,a,1x,i0)'
    character(len=30)          :: real_str

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    if (rank == 0) then
        call print_head_casci
        print '(2(A,1X,I0))', 'initialization of mpi, rank :', rank, ' nprocs :', nprocs
        print *, ''
        print *, ' ENTER R4DCASCI PROGRAM'
        print *, ''
    end if
    tmem = 0.0d+00
    if (rank == 0) then
        call write_allocated_memory_size
    end if
    call get_current_time(init_time); call print_time(init_time); start_time = init_time

    if (rank == 0) then
        print int_input_form, 'ninact        =', ninact
        print int_input_form, 'nact          =', nact
        print int_input_form, 'nsec          =', nsec
        print int_input_form, 'nelec         =', nelec
        print int_input_form, 'nroot         =', nroot
        print int_input_form, 'selectroot    =', selectroot
        print int_input_form, 'totsym        =', totsym
        write (real_str, '(E20.10)') eshift
        print '(1x,a,1x,a)', 'eshift        =', trim(adjustl(real_str))
        print int_input_form, 'diracver      =', dirac_version
        print int_input_form, 'scheme        =', mdcint_scheme
        if (ras1_size /= 0) print *, "RAS1 =", ras1_list
        if (ras2_size /= 0) print *, "RAS2 =", ras2_list
        if (ras3_size /= 0) print *, "RAS3 =", ras3_list
        print int_input_form, 'ras1_max_hole =', ras1_max_hole
        print int_input_form, 'ras3_max_elec =', ras3_max_elec
        print int_input_form, 'minholeras1   =', min_hole_ras1
        print *, 'debugprint    =', debug
        if (enable_restart) print *, "restart       =", enable_restart
        print *, ''
    end if

    ! Read around the MDCINT file and determine if the imaginary part of the 2-electron integral is written or not.
    call check_realonly()
    if (skip_mdcint) then
        if (rank == 0) print *, "Skip create_newmdcint (Activated skip_mdcint option by user input file)"
        call get_mdcint_filename(0)
    else
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
        ! Create UTChem type MDCINT file from Dirac MDCINT file
        call create_newmdcint
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
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
        print '(" IREP IS ",A6)', repna(totsym)
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

    ! Deallocate the memory
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
    if (allocated(inttwi)) then
        Call memminus(KIND(inttwi), SIZE(inttwi), 1); deallocate (inttwi)
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
    if (rank == 0) call write_allocated_memory_size
    call get_current_time_and_print_diff(init_time, end_time) ! Print the total time
    if (rank == 0) print *, 'End r4dcasci part'

END subroutine r4dcasci
