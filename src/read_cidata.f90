subroutine read_cidata
    ! Read CASCI energy and CI coefficients and dictionary of CAS configurations from cidata file
    use, intrinsic :: iso_fortran_env, only: int64
    use module_dict, only: add, destruct_dict
    use module_error, only: stop_with_errorcode
    use module_essential_input
    use module_file_manager, only: open_unformatted_file
    use module_realonly, only: realonly
    use module_global_variables
    implicit none
    character(len=cidata_key_size) :: key
    integer :: unit, i, dict_cas_idx_size
    integer :: ninact_read, nact_read, nsec_read, nelec_read, nroot_read
    integer(kind=int64), allocatable :: dict_cas_idx_values(:)
    real(8), allocatable    :: ecas(:)

    if (allocated(essential_inputs)) deallocate (essential_inputs)
    call add_essential_input("ninact")
    call add_essential_input("nact")
    call add_essential_input("nsec")
    call add_essential_input("nelec")
    call add_essential_input("ndet")
    call add_essential_input("nroot")
    call add_essential_input("ecas")
    call add_essential_input("dict_cas_idx_values")
    call add_essential_input("ci_coefficients")
    call add_essential_input("end")

    call open_unformatted_file(unit, file="CIDATA", status="old", optional_position="append")
    backspace (unit)
    read (unit) key
    if (trim(adjustl(key)) /= "end") then ! cidata must be ended with "end"
        if (rank == 0) print *, "Error: cidata file is corrupted."
        call stop_with_errorcode(1)
    end if
    rewind (unit)

    do while (.true.)
        read (unit) key
        select case (trim(adjustl(key)))
        case ("ninact")
            read (unit) ninact_read
            if (ninact_read /= ninact) then
                if (rank == 0) print *, "Error: ninact in cidata file is not equal to ninact in input file."
                call stop_with_errorcode(1)
            end if
            call update_esesential_input(trim(adjustl(key)), .true.)
        case ("nact")
            read (unit) nact_read
            if (nact_read /= nact) then
                if (rank == 0) print *, "Error: nact in cidata file is not equal to nact in input file."
                call stop_with_errorcode(1)
            end if
            call update_esesential_input(trim(adjustl(key)), .true.)
        case ("nsec")
            read (unit) nsec_read
            if (nsec_read /= nsec) then
                if (rank == 0) print *, "Error: nsec in cidata file is not equal to nsec in input file."
                call stop_with_errorcode(1)
            end if
            call update_esesential_input(trim(adjustl(key)), .true.)
        case ("nelec")
            read (unit) nelec_read
            if (nelec_read /= nelec) then
                if (rank == 0) print *, "Error: nelec in cidata file is not equal to nelec in input file."
                call stop_with_errorcode(1)
            end if
            call update_esesential_input(trim(adjustl(key)), .true.)
        case ("ndet")
            read (unit) ndet
            if (ndet < 0) then
                if (rank == 0) print *, "Error: Invalid ndet in cidata file. ndet = ", ndet
                call stop_with_errorcode(1)
            end if
            call update_esesential_input(trim(adjustl(key)), .true.)
        case ("nroot")
            read (unit) nroot_read
            if (nroot_read < 0) then
                if (rank == 0) print *, "Error: Invalid nroot in cidata file. nroot = ", nroot_read
                call stop_with_errorcode(1)
            end if
            call update_esesential_input(trim(adjustl(key)), .true.)
        case ("ecas")
            if (.not. essential_input_is_specified("nroot")) then
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
            if (.not. essential_input_is_specified("ndet")) then
                if (rank == 0) print *, "Error: dict_cas_idx_values detected before ndet."
                call stop_with_errorcode(1)
            end if
            allocate (dict_cas_idx_values(ndet))
            read (unit) dict_cas_idx_values
            call destruct_dict(dict_cas_idx)
            call destruct_dict(dict_cas_idx_reverse)
            do i = 1, ndet
                call add(dict_cas_idx, int(i, kind=int64), dict_cas_idx_values(i))
                call add(dict_cas_idx_reverse, dict_cas_idx_values(i), int(i, kind=int64))
            end do
            deallocate (dict_cas_idx_values)
            call update_esesential_input(trim(adjustl(key)), .true.)
        case ("ci_coefficients")
            if (.not. essential_input_is_specified("ndet")) then
                if (rank == 0) print *, "Error: ci_coefficients detected before ndet."
                call stop_with_errorcode(1)
            end if
            if (.not. essential_input_is_specified("nroot")) then
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
            if (rank == 0) print *, "Error: Unknown keyword in cidata file."
            call stop_with_errorcode(1)
        end select
    end do

    call check_all_essential_inputs_specified
    deallocate (essential_inputs)
    close (unit)

end subroutine read_cidata
