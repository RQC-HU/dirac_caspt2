module dcaspt2_restart_file

    use module_error, only: stop_with_errorcode
    use module_file_manager, only: open_formatted_file, check_iostat
    use module_global_variables, only: rank, enable_restart, sumc2, sumc2_subspace, e2all, e2_subspace, next_subspace, &
                                       len_convert_int_to_chr, totsym, selectroot
    implicit none

    private
    character(len=1) :: subspace_list(9) = (/'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'/)  ! I is not a valid subspace, but it is used to skip A-H subspace calc and print the final energy

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
        case ('I')
            get_subspace_idx = 9 ! I is not a valid subspace, but it is used to skip A-H subspace calc and print the final energy
        case default
            if (rank == 0) print *, "Invalid subspace (not A-H), no need to restart CASPT2 calculation"
            call stop_with_errorcode(1)
        end select
    end function get_subspace_idx

    subroutine read_and_validate_restart_file
        implicit none

        character(len=*), parameter :: key_subspace = "Subspace:", key_sumc2 = "Sumc2:", key_energy = "Energy:"
        character(len=*), parameter :: restart_file_base = "caspt2_restart"
        character(:), allocatable :: restart_file
        character(len=len_convert_int_to_chr) :: chr_totsym, chr_root
        character(len=500) :: buf, buf_internal
        character(len=1) :: subspace
        real(kind=8) :: read_sumc2, read_energy
        integer :: i, idx, file_unit, iostat, subspace_index, totsym_read, selectroot_read
        logical :: restart_file_exists, eof, valid_subspace = .false.

        write (chr_totsym, *) totsym
        write (chr_root, *) selectroot
        restart_file = trim(adjustl(restart_file_base))//"_"//trim(adjustl(chr_totsym))//"_"//trim(adjustl(chr_root))
        print *, "Restart file: ", restart_file
        inquire (file=restart_file, exist=restart_file_exists)
        if (.not. restart_file_exists) then
            if (rank == 0) print *, restart_file//" does not exist"
            return  ! No restart file, continue CASPT2 calculation from scratch
        end if

        call open_formatted_file(file_unit, restart_file, "old")
        call read_totsym(totsym_read)
        if (totsym_read /= totsym) then
            if (rank == 0) print '(2(A,I0,1x))', "totsym = ", totsym, ", restart file totsym = ", totsym_read
            call error_restart_file("totsym in "//restart_file//" does not match with the input file")
        end if
        call read_selectroot(selectroot_read)
        if (selectroot_read /= selectroot) then
            if (rank == 0) print '(2(A,I0,1x))', "selectroot = ", selectroot, ", restart file selectroot = ", selectroot_read
            call error_restart_file("selectroot in "//restart_file//" does not match with the input file")
        end if
        call read_next_subspace(subspace_index)

        do while (.true.)
            read (file_unit, '(A)', iostat=iostat) buf
            call check_iostat(iostat, restart_file, eof)
            if (eof) exit

            ! Read subspace
            idx = index(buf, key_subspace) + len(key_subspace)
            if (idx == 0) call error_restart_file("Error reading subspace info in caspt2_restart")
            read (buf(idx:), *, iostat=iostat) buf_internal
            if (iostat /= 0) call error_restart_file("Error reading subspace info in caspt2_restart")
            buf_internal = trim(adjustl(buf_internal))
            read (buf_internal, *, iostat=iostat) subspace
            if (iostat /= 0) call error_restart_file("Error reading subspace info in caspt2_restart")
            subspace_index = get_subspace_idx(subspace)

            ! Read sumc2
            idx = index(buf, key_sumc2) + len(key_sumc2)
            if (idx == 0) call error_restart_file("Error reading sumc2 info in caspt2_restart")
            read (buf(idx:), *, iostat=iostat) buf_internal
            if (iostat /= 0) call error_restart_file("Error reading sumc2 info in caspt2_restart")
            read (buf_internal, *, iostat=iostat) read_sumc2
            if (iostat /= 0) call error_restart_file("Error reading sumc2 info in caspt2_restart")
            sumc2_subspace(subspace_index) = read_sumc2
            sumc2 = sumc2 + read_sumc2

            ! Read energy
            idx = index(buf, key_energy) + len(key_energy)
            if (idx == 0) call error_restart_file("Error reading energy info in caspt2_restart")
            read (buf(idx:), *, iostat=iostat) buf_internal
            if (iostat /= 0) call error_restart_file("Error reading energy info in caspt2_restart")
            read (buf_internal, *, iostat=iostat) read_energy
            if (iostat /= 0) call error_restart_file("Error reading energy info in caspt2_restart")
            e2_subspace(subspace_index) = read_energy
            e2all = e2all + read_energy
        end do

    contains
        subroutine read_totsym(totsym_ret)
            implicit none
            integer, intent(out) :: totsym_ret
            read (file_unit, '(A)', iostat=iostat) buf
            call check_iostat(iostat, restart_file, eof)
            idx = index(buf, ":")
            if (idx == 0) call error_restart_file("Error reading totsym info in caspt2_restart")
            read (buf(idx + 1:), '(A)', iostat=iostat) buf_internal
            if (iostat /= 0) call error_restart_file("Error reading totsym info in caspt2_restart")
            read (buf_internal, *) totsym_ret ! convert string to integer
        end subroutine read_totsym

        subroutine read_selectroot(selectroot_ret)
            implicit none
            integer, intent(out) :: selectroot_ret
            read (file_unit, '(A)', iostat=iostat) buf
            call check_iostat(iostat, restart_file, eof)
            idx = index(buf, ":")
            if (idx == 0) call error_restart_file("Error reading selectroot info in caspt2_restart")
            read (buf(idx + 1:), '(A)', iostat=iostat) buf_internal
            if (iostat /= 0) call error_restart_file("Error reading selectroot info in caspt2_restart")
            read (buf_internal, *) selectroot_ret ! convert string to integer
        end subroutine read_selectroot

        subroutine read_next_subspace(subpsace_idx_ret)
            implicit none
            integer, intent(out) :: subpsace_idx_ret
            read (file_unit, '(A)', iostat=iostat) buf
            call check_iostat(iostat, restart_file, eof)
            idx = index(buf, ":")
            if (idx == 0) call error_restart_file("Error reading next subspace info in caspt2_restart")
            read (buf(idx + 1:), '(A)', iostat=iostat) buf_internal
            if (iostat /= 0) call error_restart_file("Error reading next subspace info in caspt2_restart")

            buf_internal = trim(adjustl(buf_internal))
            read (buf_internal, *) next_subspace
            subpsace_idx_ret = get_subspace_idx(next_subspace)
        end subroutine read_next_subspace

        subroutine error_restart_file(error_message)
            character(len=*), intent(in) :: error_message
            if (rank == 0) print *, error_message
            call stop_with_errorcode(1)
        end subroutine error_restart_file
    end subroutine read_and_validate_restart_file
end module dcaspt2_restart_file
