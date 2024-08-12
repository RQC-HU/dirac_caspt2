module dcaspt2_restart_file

    use module_error, only: stop_with_errorcode
    use module_file_manager, only: open_formatted_file, check_iostat
    use module_global_variables, only: rank, enable_restart, sumc2, sumc2_subspace, e2all, e2_subspace, next_subspace
    implicit none

    private
    character(len=1) :: subspace_list(8) = (/'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'/)

    public :: is_skip_restart_file_subspace, get_subspace_idx, read_and_validate_restart_file
contains

    function is_skip_restart_file_subspace(subspace_char)
        character(len=1), intent(in) :: subspace_char
        logical :: is_skip_restart_file_subspace
        is_skip_restart_file_subspace = .FALSE.
        if (.not. enable_restart) return  ! Don't skip

        if (get_subspace_idx(next_subspace) > get_subspace_idx(subspace_char)) then
            is_skip_restart_file_subspace = .TRUE.
        end if
    end function is_skip_restart_file_subspace

    function get_subspace_idx(subspace_char)
        character(len=1), intent(in) :: subspace_char
        integer :: get_subspace_idx

        select case (subspace_char)
        case ('A')
            get_subspace_idx = 1
        case ('B')
            get_subspace_idx = 2
        case ('C')
            get_subspace_idx = 3
        case ('D')
            get_subspace_idx = 4
        case ('E')
            get_subspace_idx = 5
        case ('F')
            get_subspace_idx = 6
        case ('G')
            get_subspace_idx = 7
        case ('H')
            get_subspace_idx = 8
        case default
            if (rank == 0) print *, "Invalid subspace (not A-H), no need to restart CASPT2 calculation"
            call stop_with_errorcode(1)
        end select
    end function get_subspace_idx

    subroutine read_and_validate_restart_file
        implicit none

        character(len=*), parameter :: key_subspace = "Subspace:", key_sumc2 = "Sumc2:", key_energy = "Energy:"
        character(len=*), parameter :: restart_file = "caspt2_restart"
        character(len=500) :: buf, buf_internal
        character(len=1) :: subspace
        real(kind=8) :: read_sumc2, read_energy
        integer :: i, idx, file_unit, iostat, subspace_index
        logical :: restart_file_exists, eof, valid_subspace = .false.

        inquire (file=restart_file, exist=restart_file_exists)
        if (.not. restart_file_exists) call error_restart_file("caspt2_restart file does not exist")

        call open_formatted_file(file_unit, "caspt2_restart", "old")
        read (file_unit, '(A)', iostat=iostat) buf
        call check_iostat(iostat, restart_file, eof)
        idx = index(buf, ":")
        if (idx == 0) call error_restart_file("Error reading next subspace info in caspt2_restart")
        read (buf(idx + 1:), '(A)', iostat=iostat) buf_internal
        if (iostat /= 0) call error_restart_file("Error reading next subspace info in caspt2_restart")
        buf_internal = trim(adjustl(buf_internal))
        read (buf_internal, *) next_subspace
        subspace_index = get_subspace_idx(next_subspace)

        do while (.true.)
            read (file_unit, '(A)', iostat=iostat) buf
            call check_iostat(iostat, restart_file, eof)
            if (eof) exit

            ! Read subspace
            idx = index(buf, key_subspace) + len(key_subspace)
            read (buf(idx:), *) buf_internal
            buf_internal = trim(adjustl(buf_internal))
            read (buf_internal, *, iostat=iostat) subspace
            if (iostat /= 0) call error_restart_file("Error reading subspace info in caspt2_restart")
            subspace_index = get_subspace_idx(subspace)

            ! Read sumc2
            idx = index(buf, key_sumc2) + len(key_sumc2)
            read (buf(idx:), *, iostat=iostat) buf_internal
            if (iostat /= 0) call error_restart_file("Error reading sumc2 info in caspt2_restart")
            read (buf_internal, *, iostat=iostat) read_sumc2
            if (iostat /= 0) call error_restart_file("Error reading sumc2 info in caspt2_restart")
            sumc2_subspace(subspace_index) = read_sumc2
            sumc2 = sumc2 + read_sumc2

            ! Read energy
            idx = index(buf, key_energy) + len(key_energy)
            read (buf(idx:), *, iostat=iostat) buf_internal
            if (iostat /= 0) call error_restart_file("Error reading energy info in caspt2_restart")
            read (buf_internal, *, iostat=iostat) read_energy
            if (iostat /= 0) call error_restart_file("Error reading energy info in caspt2_restart")
            e2_subspace(subspace_index) = read_energy
            e2all = e2all + read_energy
        end do

    contains
        subroutine error_restart_file(error_message)
            character(len=*), intent(in) :: error_message
            if (rank == 0) print *, error_message
            call stop_with_errorcode(1)
        end subroutine error_restart_file
    end subroutine read_and_validate_restart_file
end module dcaspt2_restart_file
