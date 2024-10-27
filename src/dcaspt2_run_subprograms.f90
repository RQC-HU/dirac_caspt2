subroutine dcaspt2_run_subprograms
    use dcaspt2_restart_file, only: read_and_validate_restart_file
    use module_file_manager, only: open_formatted_file
    use module_global_variables
    use read_input_module, only: read_input
    implicit none
    integer :: unit_input
    character(:), allocatable :: filename

    if (rank == 0) print *, "Starting DCASPT2 calculation"

    call open_formatted_file(unit=unit_input, file='active.inp', status="old", optional_action="read")
    rewind (unit_input)
    call read_input(unit_input)
    if (enable_restart) call read_and_validate_restart_file

    ! Read MRCONEE file (orbital energies, symmetries and multiplication tables)
    if (rank == 0) print *, ' '
    if (rank == 0) print *, 'Reading MRCONEE (1-e integrals)'
    allocate(filename, source='MRCONEE')
    call check_dirac_integer_size(filename)
    call read_mrconee(filename)

    if (doivo) then
        call r4divo_co
    end if

    if (docasci) then
        call r4dcasci
    end if

    if (docaspt2) then
        call r4dcaspt2_tra
    end if

    ! Deallocate memory
    if (allocated(space_idx)) then
        Call memminus(KIND(space_idx), SIZE(space_idx), 1); deallocate (space_idx)
    end if
    if (allocated(MULTB_S)) then
        Call memminus(KIND(MULTB_S), SIZE(MULTB_S), 1); deallocate (MULTB_S)
    end if
    if (allocated(MULTB_D)) then
        Call memminus(KIND(MULTB_D), SIZE(MULTB_D), 1); deallocate (MULTB_D)
    end if
    if (allocated(MULTB_DS)) then
        Call memminus(KIND(MULTB_DS), SIZE(MULTB_DS), 1); deallocate (MULTB_DS)
    end if
    if (allocated(irpamo)) then
        Call memminus(KIND(irpamo), SIZE(irpamo), 1); deallocate (irpamo)
    end if
    if (allocated(indmo_cas_to_dirac)) then
        Call memminus(KIND(indmo_cas_to_dirac), SIZE(indmo_cas_to_dirac), 1); deallocate (indmo_cas_to_dirac)
    end if
    if (allocated(indmo_dirac_to_cas)) then
        Call memminus(KIND(indmo_dirac_to_cas), SIZE(indmo_dirac_to_cas), 1); deallocate (indmo_dirac_to_cas)
    end if
    if (allocated(one_elec_int_i)) then
        Call memminus(KIND(one_elec_int_i), SIZE(one_elec_int_i), 1); deallocate (one_elec_int_i)
    end if
    if (allocated(one_elec_int_r)) then
        Call memminus(KIND(one_elec_int_r), SIZE(one_elec_int_r), 1); deallocate (one_elec_int_r)
    end if

end subroutine dcaspt2_run_subprograms
