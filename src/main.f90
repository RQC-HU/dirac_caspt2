program main

    use module_file_manager, only: open_formatted_file
    use module_global_variables
    use module_mpi
    use read_input_module, only: read_input
    use module_realonly, only: check_realonly
    implicit none
    integer :: unit_input, idx_totsym, idx_selectroot

    call mpi_init_wrapper
    if (rank == 0) then
        print '(2(A,1X,I0))', 'initialization of mpi, rank :', rank, ' nprocs :', nprocs
        print *, ''
        print *, ' ENTER DIRAC-CASPT2 PROGRAM written by M. Abe 2007.7.19'
        print *, ''
    end if

    tmem = 0.0d+00

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
    if (rank == 0) print *, ' ENTER READ MRCONEE'
    call read_mrconee('MRCONEE')
    ! Read around the MDCINT file and determine if the imaginary part of the 2-electron integral is written or not.
    call check_realonly()

    print *, "totsym_list, selectroot_list", totsym_list, selectroot_list
    do idx_totsym = 1, size(totsym_list)
        totsym = totsym_list(idx_totsym)
        do idx_selectroot = 1, size(selectroot_list, 1)
            selectroot = selectroot_list(idx_selectroot, idx_totsym)
            if (selectroot == 0) exit
            call r4dcasci
            call mpi_barrier_wrapper
            call r4dcaspt2_tra
            call mpi_barrier_wrapper
            print *, "totsym, selectroot", totsym, selectroot, idx_selectroot, idx_totsym
        end do
        call mpi_barrier_wrapper
    end do
    call mpi_finalize_wrapper

contains
    subroutine deallocate_variables
        implicit none

    end subroutine deallocate_variables

end program main
