subroutine dcaspt2_run_subprograms
    use dcaspt2_restart_file, only: read_and_validate_restart_file
    use module_file_manager, only: open_formatted_file
    use module_global_variables
    use module_realonly, only: check_realonly
    use read_input_module, only: read_input
    implicit none
    integer :: unit_input, i
    character(:), allocatable :: filename

    call print_head

    call open_formatted_file(unit=unit_input, file='active.inp', status="old", optional_action="read")
    rewind (unit_input)
    call read_input(unit_input)

    ! Read MRCONEE file (orbital energies, symmetries and multiplication tables)
    if (rank == 0) then
        print *, ' '
        print *, 'Reading MRCONEE (1-e integrals)'
    end if
    allocate (filename, source='MRCONEE')
    call check_dirac_integer_size(filename)
    call read_mrconee(filename)
    call check_realonly

    if (docountndet) then
        call search_cas_configuration ! If docountndet is true, only search the number of CASCI configurations and print them
        return ! skip subprograms if docountndet is true
    end if

    if (doivo) then
        call r4divo_co
        call dcaspt2_deallocate
        return
    end if

    do i = 1, size(caspt2_ciroots, 1)
        totsym = caspt2_ciroots(i, 1)
        selectroot = caspt2_ciroots(i, 2)
        nroot = max_selectroot_list(totsym)
        e2_subspace = 0
        sumc2 = 0
        sumc2_subspace = 0
        if (enable_restart) call read_and_validate_restart_file

        if (docasci .and. .not. casci_done(totsym)) then
            call r4dcasci
            casci_done(totsym) = .true.
        end if

        if (docaspt2) then
            call r4dcaspt2_tra
        end if
    end do
    call dcaspt2_deallocate

end subroutine dcaspt2_run_subprograms
