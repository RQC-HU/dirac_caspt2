subroutine read_cimat
    ! Read CASCI energy and CI coefficients and dictionary of CAS configurations from CIMAT file
    use module_dict, only: add
    use module_error, only: stop_with_errorcode
    use module_essential_input
    use module_file_manager, only: open_unformatted_file
    use module_realonly, only: realonly
    use module_global_variables
    implicit none
    character(len=cimat_key_size) :: key
    integer :: unit, i, dict_cas_idx_size, nroot_read
    integer(8), allocatable :: dict_cas_idx_values(:)
    real(8), allocatable    :: ecas(:)

    ! Initialize variables as invalid values
    dict_cas_idx_size = -1
    ndet = -1
    nroot_read = -1

    if (allocated(essential_inputs)) deallocate (essential_inputs)
    call add_essential_input("ndet")
    call add_essential_input("nroot")
    call add_essential_input("ecas")
    call add_essential_input("dict_cas_idx_values")
    call add_essential_input("ci_coefficients")
    call add_essential_input("end")

    call open_unformatted_file(unit, "CIMAT", "old", "read", "append")
    backspace (unit)
    read (unit) key
    if (trim(adjustl(key)) /= "end") then
        if (rank == 0) print *, "Error: CIMAT file is corrupted."
        call stop_with_errorcode(1)
    end if
    rewind (unit)

    do while (.true.)
        read (unit) key
        select case (trim(adjustl(key)))
        case ("ndet")
            read (unit) ndet
            if (ndet < 0) then
                if (rank == 0) print *, "Error: Invalid ndet in CIMAT file. ndet = ", ndet
                call stop_with_errorcode(1)
            end if
            call update_esesential_input(trim(adjustl(key)), .true.)
        case ("nroot")
            read (unit) nroot_read
            if (nroot_read < 0) then
                if (rank == 0) print *, "Error: Invalid nroot in CIMAT file. nroot = ", nroot_read
                call stop_with_errorcode(1)
            end if
            call update_esesential_input(trim(adjustl(key)), .true.)
        case ("ecas")
            if (nroot_read == -1) then
                if (rank == 0) print *, "Error: ecas detected before nroot."
                call stop_with_errorcode(1)
            end if
            allocate (ecas(nroot_read))
            allocate (eigen(nroot_read)); Call memplus(KIND(eigen), SIZE(eigen), 1)
            ! Read CASCI energy
            read (unit) ecas
            eigen(:) = ecas(1:nroot_read) + ecore
            deallocate (ecas)
            call update_esesential_input(trim(adjustl(key)), .true.)
        case ("dict_cas_idx_values")
            if (ndet == -1) then
                if (rank == 0) print *, "Error: dict_cas_idx_values detected before ndet."
                call stop_with_errorcode(1)
            end if
            allocate (dict_cas_idx_values(ndet))
            read (unit) dict_cas_idx_values
            do i = 1, ndet
                call add(dict_cas_idx, i, dict_cas_idx_values(i))
                call add(dict_cas_idx_reverse, dict_cas_idx_values(i), i)
            end do
            call update_esesential_input(trim(adjustl(key)), .true.)
        case ("ci_coefficients")
            if (ndet == -1) then
                if (rank == 0) print *, "Error: ci_coefficients detected before ndet."
                call stop_with_errorcode(1)
            end if
            if (nroot_read == -1) then
                if (rank == 0) print *, "Error: ci_coefficients detected before nroot."
                call stop_with_errorcode(1)
            end if
            allocate (cir(ndet, nroot_read), cii(ndet, nroot_read))
            call memplus(KIND(cir), SIZE(cir), 1); call memplus(KIND(cii), SIZE(cii), 1)
            if (realonly%is_realonly()) then
                read (unit) cir
                cii = 0.0d0
            else
                read (unit) cir
                read (unit) cii
            end if
            call update_esesential_input(trim(adjustl(key)), .true.)
        case ("end")
            call update_esesential_input(trim(adjustl(key)), .true.)
            exit
        case default
            if (rank == 0) print *, "Error: Unknown keyword in CIMAT file."
            call stop_with_errorcode(1)
        end select
    end do

    call check_all_essential_inputs_specified
    deallocate (essential_inputs)

end subroutine read_cimat
