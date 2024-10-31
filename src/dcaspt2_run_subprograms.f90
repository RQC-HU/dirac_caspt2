subroutine dcaspt2_run_subprograms
    use dcaspt2_restart_file, only: read_and_validate_restart_file
    use module_file_manager, only: open_formatted_file
    use module_global_variables
    use module_realonly, only: check_realonly
    use read_input_module, only: read_input
    implicit none
    integer :: unit_input
    character(:), allocatable :: filename

    call print_head

    call open_formatted_file(unit=unit_input, file='active.inp', status="old", optional_action="read")
    rewind (unit_input)
    call read_input(unit_input)
    if (enable_restart) call read_and_validate_restart_file

    ! Read MRCONEE file (orbital energies, symmetries and multiplication tables)
    if (rank == 0) print *, ' '
    if (rank == 0) print *, 'Reading MRCONEE (1-e integrals)'
    allocate (filename, source='MRCONEE')
    call check_dirac_integer_size(filename)
    call read_mrconee(filename)
    call check_realonly

    if (doivo) then
        call r4divo_co
    end if

    if (docasci) then
        call r4dcasci
    end if

    if (docaspt2) then
        call r4dcaspt2_tra
    end if

    call dcaspt2_deallocate

end subroutine dcaspt2_run_subprograms
