module read_input_module
!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
! read_input_module
! Copyright (c) by the authors of DIRAC-CASPT2.
! Author K.Noda
!
! This is a utility module that interpret and parse input strings.
!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
    use module_essential_input
    use module_global_variables, only: rank, len_convert_int_to_chr
    use module_error, only: stop_with_errorcode
    implicit none
    private
    public read_input, check_substring, ras_read, lowercase, uppercase
    logical is_end, set_caspt2_ciroots
    integer, parameter :: input_intmax = 10**9, max_str_length = 500
    type(essential_inputs_container) :: container

contains

    subroutine init_essential_variables
        call container%add_essential_input(".ninact")
        call container%add_essential_input(".nact")
        call container%add_essential_input(".nsec")
        call container%add_essential_input(".nelec")
        call container%add_essential_input(".diracver")
        call container%add_essential_input(".subprograms")
    end subroutine init_essential_variables

    subroutine print_input_file(unit_num)
        implicit none
        integer, intent(in) :: unit_num
        character(len=max_str_length) :: line
        integer :: iostat

        if (rank == 0) then
            print *, ""
            print *, "Input file:"
            print *, "```inp"
            do while (.true.)
                read (unit_num, '(A)', iostat=iostat) line
                if (iostat /= 0) exit
                print *, trim(adjustl(line))
            end do
            print *, "```"
            print *, ""
        end if
    end subroutine print_input_file

    subroutine read_input(unit_num, bypass_reqired_file_check)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine is the entry point to read active.inp
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        use module_global_variables, only: ras1_size, ras2_size, ras3_size
        use module_index_utils, only: set_global_index
        implicit none
        integer, intent(in) :: unit_num
        logical, optional, intent(in) :: bypass_reqired_file_check ! It is used for unit test
        integer :: iostat
        character(len=max_str_length) :: string
        logical :: is_comment, do_reqired_file_check
        is_end = .false.
        set_caspt2_ciroots = .false.

        if (present(bypass_reqired_file_check)) then
            do_reqired_file_check = .not. bypass_reqired_file_check
        else
            do_reqired_file_check = .true.
        end if

        call print_input_file(unit_num)
        rewind (unit_num)
        call init_essential_variables

        do while (.not. is_end) ! Read the input file until the "end" is found
            read (unit_num, "(a)", iostat=iostat) string
            if (iostat < 0) then
                if (rank == 0) print *, "ERROR: YOU NEED TO ADD 'end' AT THE END OF YOUR INPUT FILE."
                call stop_with_errorcode(1)
            else if (iostat > 0) then
                if (rank == 0) print *, "ERROR: Error in input, failed to read active.inp"
                call stop_with_errorcode(1)
            end if
            call is_comment_line(string, is_comment)
            if (is_comment) cycle ! Read the next line
            call read_keyword_and_value(unit_num, string)
        end do

        call container%check_all_essential_inputs_specified()
        call set_global_index
        call set_mdcint_scheme
        call check_ciroots_set
        if (do_reqired_file_check) call check_reqired_files_exist
        ! Check the RAS configuration
        if (ras1_size /= 0 .or. ras2_size /= 0 .or. ras3_size /= 0) call check_ras_is_valid

    end subroutine read_input

    subroutine read_keyword_and_value(unit_num, string)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine recognize the type of input that follows from the next line
        ! and calls the subroutine that we must call
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        use module_global_variables
        implicit none
        integer, intent(in) :: unit_num
        character(*), intent(inout) :: string
        character(len=max_str_length) :: input
        logical :: is_comment
        call lowercase(string)
        select case (trim(adjustl(string)))

        case (".ninact")
            call read_an_integer(unit_num, ".ninact", 0, input_intmax, ninact)
            call container%update_essential_input(".ninact", .true.)

        case (".nact")
            call read_an_integer(unit_num, ".nact", 0, input_intmax, nact)
            call container%update_essential_input(".nact", .true.)

        case (".nsec")
            call read_an_integer(unit_num, ".nsec", 0, input_intmax, nsec)
            call container%update_essential_input(".nsec", .true.)

        case (".nelec")
            call read_an_integer(unit_num, ".nelec", 0, input_intmax, nelec)
            call container%update_essential_input(".nelec", .true.)

        case (".caspt2_ciroots")
            call read_caspt2_ciroots(unit_num)
            set_caspt2_ciroots = .true.

        case (".ncore")
            call read_an_integer(unit_num, ".ncore", 0, input_intmax, ncore)

        case (".eshift")
            do
                read (unit_num, '(A)') input
                call is_comment_line(input, is_comment)
                if (.not. is_comment) then
                    read (input, *) eshift
                    exit
                end if
            end do

        case (".diracver")
            call read_an_integer(unit_num, ".diracver", 0, input_intmax, dirac_version)
            call container%update_essential_input(".diracver", .true.)

        case (".nhomo")
            call read_an_integer(unit_num, ".nhomo", 0, input_intmax, nhomo)

        case (".ras1")
            call ras_read(unit_num, ras1_list, 1)
            ras1_size = size(ras1_list, 1)
            call read_an_integer(unit_num, ".ras1", 0, ras1_size, ras1_max_hole)

        case (".ras2")
            call ras_read(unit_num, ras2_list, 2)
            ras2_size = size(ras2_list, 1)

        case (".ras3")
            call ras_read(unit_num, ras3_list, 3)
            ras3_size = size(ras3_list, 1)
            call read_an_integer(unit_num, ".ras3", 0, ras3_size, ras3_max_elec)

        case (".minholeras1")
            call read_an_integer(unit_num, ".minholeras1", 0, input_intmax, min_hole_ras1)

        case (".minelecras3")
            call read_an_integer(unit_num, ".minelecras3", 0, input_intmax, min_elec_ras3)

        case (".nocc")
            if (inversion) call err_ivo_input
            call read_an_integer(unit_num, ".nocc", 0, input_intmax, occ_mo_num(1))
            no_inversion = .true.

        case (".noccg")
            if (no_inversion) call err_ivo_input
            call read_an_integer(unit_num, ".noccg", 0, input_intmax, occ_mo_num(1))
            inversion = .true.

        case (".noccu")
            if (no_inversion) call err_ivo_input
            call read_an_integer(unit_num, ".noccu", 0, input_intmax, occ_mo_num(2))
            inversion = .true.

        case (".nvcut")
            if (inversion) call err_ivo_input
            call read_an_integer(unit_num, ".nvcut", 0, input_intmax, vcut_mo_num(1))
            no_inversion = .true.

        case (".nvcutg")
            if (no_inversion) call err_ivo_input
            call read_an_integer(unit_num, ".nvcutg", 0, input_intmax, vcut_mo_num(1))
            inversion = .true.

        case (".nvcutu")
            if (no_inversion) call err_ivo_input
            call read_an_integer(unit_num, ".nvcutu", 0, input_intmax, vcut_mo_num(2))
            inversion = .true.

        case (".scheme")
            call read_an_integer(unit_num, ".scheme", 1, input_intmax, mdcint_scheme)
            is_scheme_set = .true.

        case (".debugprint")
            debug = .true.

        case (".restart")
            enable_restart = .true.

        case (".subprograms")
            call read_subprograms(unit_num)
            call container%update_essential_input(".subprograms", .true.)

        case (".countndet")
            docountndet = .true.
            ! .countndet just calls the search_cas_configuration subroutine, so it is not a subprogram,
            ! but if .countndet is specified, the other subroutines will be skipped.
            ! Therefore, if .countndet is specified, .subprograms doesn't need to be specified.
            ! Thus, we set essential input "subprograms" to .true.
            call container%update_essential_input(".subprograms", .true.)

        case (".end")
            is_end = .true.

        case default
            if (rank == 0) print *, "ERROR: Unknown input: ", trim(string)
            call stop_with_errorcode(1)
        end select
    contains
        subroutine err_ivo_input
            implicit none
            if (rank == 0) print *, "ERROR: nocc or nvcut and noccg or noccu or nvcutg or nvcutu", &
                "cannot be specified at the same time."
            call stop_with_errorcode(1)
        end subroutine err_ivo_input

    end subroutine read_keyword_and_value

    subroutine ras_read(unit_num, ras_list, ras_num)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine returns RAS[1,2,3] list from the user input
        ! (e.g.) INPUT  : string = "1,2,4..10,13,17..20"
        !        OUTPUT : ras[1,2,3]_list = [1,2,4,5,6,7,8,9,10,13,17,18,19,20], (ras[1,2,3]_list is a global list in module_global_variables)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        use module_global_variables, only: max_ras_spinor_num
        use module_sort_swap, only: heapSort
        implicit none
        integer, allocatable, intent(inout) :: ras_list(:)
        integer, intent(in) :: unit_num, ras_num
        character(len=max_str_length) :: tmp_ras_chr
        character(:), allocatable :: ras_chr
        character(len=max_str_length) :: string
        integer :: tmp_ras(max_ras_spinor_num), idx_filled, iostat, idx

        ! store ras_num character to ras_chr
        write (tmp_ras_chr, *) ras_num
        allocate (ras_chr, source=trim(adjustl(tmp_ras_chr)))

        read (unit_num, '(a)', iostat=iostat) string ! Read a line of active.inp
        if (iostat /= 0) then
            if (rank == 0) print *, "ERROR: ras_read: iostat = ", iostat, ", string =", string
            call stop_with_errorcode(iostat)
            call exit(iostat)
        end if
        idx_filled = 0

        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        !  Parse the expressions of the form a..b in the input string and expands it to a list.
        !  (e.g.) INPUT  : string = "1,3,5..8,10", tmp_ras = [0,0,...,0],                idx_filled = 0
        !         OUTPUT : string = " , ,    ,  ", tmp_ras = [5,6,7,8,1,3,10,0,0,...,0], idx_filled = 7
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        call parse_input_string_to_int_list(string=string, list=tmp_ras, filled_num=idx_filled, &
                                            allow_int_min=0, allow_int_max=input_intmax)
        ! Does the input string contain at least one varible?
        if (idx_filled <= 0) then
            print *, "ERROR: string:", string, " rank:", rank
            call write_error_and_stop_ras_read
        end if
        allocate (ras_list(idx_filled))
        ras_list(:) = tmp_ras(1:idx_filled)
        call heapSort(list=ras_list, is_descending_order=.false.) ! Sort the ras_list in ascending order (lower to higher)

        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! Check the specification of input is kramers pair?
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!

        ! The size of ras_list must be even.
        if (mod(size(ras_list), 2) /= 0) then
            if (rank == 0) print *, "ERROR: The number of ras_list is not even."
            call write_error_and_stop_ras_read
        end if

        ! ras_list(idx) (idx : odd) must be odd number and equal to ras_list(idx+1) (idx : even)
        do idx = 1, size(ras_list, 1), 2
            ! Check the ras_list(idx) (idx : odd)  is odd number?
            if (mod(ras_list(idx), 2) /= 1) then
                if (rank == 0) then
                    print *, "ERROR: sorted ras", ras_chr, "(idx) (idx : odd) must be odd number."
                    print *, "idx,ras", ras_chr, "(idx) :", idx, ras_list(idx)
                    print *, "ras", ras_chr, "list :", ras_list
                end if
                call write_error_and_stop_ras_read
            end if
            ! Check the ras_list(idx+1) (idx : even) is equal to ras_list(idx) + 1 (idx : odd)?
            if (ras_list(idx) + 1 /= ras_list(idx + 1)) then
                if (rank == 0) print *, "ERROR: The ras", ras_chr, " is not kramers pair. idx,ras", ras_chr, "(idx),", &
                    "ras", ras_chr, "(idx+1) :", idx, ras_list(idx), ras_list(idx + 1)
                call write_error_and_stop_ras_read
            end if
        end do
    contains
        subroutine write_error_and_stop_ras_read
            implicit none
            print *, "ERROR: Error in input, can't read ras", ras_chr, " value!!. Stop the program. rank:", rank
            call stop_with_errorcode(1)
        end subroutine write_error_and_stop_ras_read
    end subroutine ras_read

    subroutine parse_input_string_to_int_list(string, list, filled_num, allow_int_min, allow_int_max)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine returns a list of integers
        ! It finds expressions of integer or the form a..b in the input string and expands it to a list.
        ! (e.g.) INPUT  : string = "1,2,4..10,13,17..20", list = [0,0,...,0],                                   filled_num = 0
        !        OUTPUT : string = " , ,     ,  ,      ", list = [1,2,4,5,6,7,8,9,10,13,17,18,19,20,0,0,...,0], filled_num = 14
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        implicit none
        integer, intent(in) :: allow_int_min, allow_int_max
        character(*), intent(inout) :: string
        integer, intent(inout) :: filled_num ! The number of numbers already filled in list
        integer, intent(inout) :: list(:)

        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! First, we call this subroutine to detect the expression of the form a..b and expand it to a list.
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        call parse_range_input_int(string, list, filled_num, allow_int_min, allow_int_max)

        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! Now, input string is seems to contain only space or , or integer
        ! So we parse integer from input string and expand it to a list
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        call parse_input_int(string, list, filled_num, allow_int_min, allow_int_max)

    end subroutine parse_input_string_to_int_list

    subroutine parse_input_int(string, list, filled_num, allow_int_min, allow_int_max)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine returns a list of integers
        ! It finds expressions of integer in the input string and expands it to a list.
        ! This subroutine can't detect expressions of the form a..b,
        ! so You "must" call this subroutine after call "parse_range_input_int" subroutine.
        ! (e.g.) INPUT  : string = "1,2,     ,13,      ", list = [4,5,6,7,8,9,10,17,18,19,20,0,...,0],        filled_num = 11
        !        OUTPUT : string = " , ,     ,  ,      ", list = [4,5,6,7,8,9,10,17,18,19,20,1,2,13,0,...,0], filled_num = 14
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        implicit none
        integer, intent(in) :: allow_int_min, allow_int_max
        character(*), intent(inout) :: string
        integer, intent(inout) :: filled_num ! The number of numbers already filled in list
        integer, intent(inout) :: list(:)
        character(len_convert_int_to_chr) :: min_str, max_str, read_int_str
        character(:), allocatable  :: pattern, invalid_input_message, subroutine_name
        integer :: read_int, read_int_digit, idx, iostat
        logical :: is_valid

        allocate (subroutine_name, source="parse_input_int")
        call check_args

        ! The variable idx is the first index that does not contain space or ,
        idx = verify(string, ' ,') ! idx is 0 if all characters in string are space or ,
        do while (idx /= 0)
            call check_valid_input

            ! Now we can get the num (e.g. "12,15" -> 12)
            call read_integer(read_int)
            call check_read_int_is_valid(read_int)
            call fill_number_to_list(read_int)

            ! Rewrite the section we read as blank. (e.g.) "10,3,5" -> "  ,3,5"
            write (read_int_str, *) read_int
            read_int_digit = len(trim(adjustl(read_int_str)))
            string(idx:idx + read_int_digit - 1) = ""

            idx = verify(string, ' ,')  ! next idx
        end do
    contains
        subroutine check_args
            implicit none
            call check_range_allow_int(allow_int_min, allow_int_max, is_valid)
            if (.not. is_valid) call write_parse_error_and_stop(subroutine_name, string)
            call create_valid_pattern(allow_int_min, allow_int_max, pattern, invalid_input_message)
        end subroutine check_args

        subroutine check_valid_input
            implicit none
            call check_substring(string(idx:idx), pattern, is_valid)
            if (.not. is_valid) then ! The number is NOT a integer or invalid input.
                if (rank == 0) print *, invalid_input_message, string(idx:)
                call write_parse_error_and_stop(subroutine_name, string)
            end if
        end subroutine check_valid_input

        subroutine read_integer(read_value)
            implicit none
            integer, intent(out) :: read_value
            read (string(idx:), *, iostat=iostat) read_value
            if (iostat /= 0) then
                if (rank == 0) then
                    print *, "Error in the section in reading the number, iostat = ", iostat, &
                        ", string = ", string(idx:)
                end if
                call write_parse_error_and_stop(subroutine_name, string)
            end if
        end subroutine read_integer

        subroutine check_read_int_is_valid(read_value)
            implicit none
            integer, intent(in) :: read_value
            call is_in_range_int(read_value, allow_int_min, allow_int_max, is_valid)
            if (.not. is_valid) then ! read_value is out of range
                write (min_str, *) allow_int_min
                write (max_str, *) allow_int_max
                if (rank == 0) print *, "ERROR: read_int is out of range,", &
                    "[", trim(adjustl(min_str)), ",", trim(adjustl(max_str)), "]"
                call write_parse_error_and_stop(subroutine_name, string)
            end if
        end subroutine check_read_int_is_valid

        subroutine fill_number_to_list(read_value)
            implicit none
            integer, intent(in) :: read_value
            if (size(list, 1) < filled_num + 1) then
                ! Can't fill numbers because the size of the list
                if (rank == 0) print *, "Can't fill the number because of the size of the list. size:", size(list, 1)
                call write_parse_error_and_stop(subroutine_name, string)
            end if
            ! Now we can store the number
            filled_num = filled_num + 1
            list(filled_num) = read_value
        end subroutine fill_number_to_list
    end subroutine parse_input_int

    subroutine parse_range_input_int(string, list, filled_num, allow_int_min, allow_int_max)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine returns a list of integers
        ! It finds expressions of the form a..b in the input string and expands it to a list.
        ! (e.g.) INPUT  : string = "1,2,4..10,13,17..20", list = [0,0,...,0],                            filled_num = 0
        !        OUTPUT : string = "1,2,     ,13,      ", list = [4,5,6,7,8,9,10,17,18,19,20,0,0,...,0], filled_num = 11
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        implicit none
        integer, intent(in) :: allow_int_min, allow_int_max
        character(*), intent(inout) :: string
        integer, intent(inout) :: filled_num ! The number of numbers already filled in list
        integer, intent(inout) :: list(:)
        character(len=len_convert_int_to_chr)     :: right_str, min_str, max_str
        character(:), allocatable   :: pattern, invalid_input_message, subroutine_name
        integer :: first_dot_index, stat, rightnum_first_idx, leftnum_idx, leftnum, rightnum, rightnum_digit, idx, iostat
        logical :: is_valid

        allocate (subroutine_name, source="parse_range_input_int")
        call check_args

        first_dot_index = index(string, '..')

        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! Read the numbers of the left and right and fill in the list to the range of the readings
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        do while (first_dot_index /= 0)

            call search_rightnum_idx(rightnum_first_idx)
            call read_rightnum(rightnum) ! read the right num (e.g. "12..15" -> 15)

            call search_leftnum_idx(leftnum_idx)
            call read_leftnum(leftnum) ! read the left num (e.g. "12..15" -> 12)

            call fill_numbers_to_list(leftnum, rightnum)

            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! Rewrite the section we read as blank. (e.g. "1, 2, 4..10, 13" -> "1, 2,      , 13" )
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            write (right_str, *) rightnum
            rightnum_digit = len(trim(adjustl(right_str))) ! Get the digit of rightnum (e.g. -10 -> 3, 23 -> 2)
            string(leftnum_idx:rightnum_first_idx + rightnum_digit - 1) = ""

            first_dot_index = index(string, '..') ! Next double dots index in string
        end do

    contains
        subroutine check_args
            implicit none
            call check_range_allow_int(allow_int_min, allow_int_max, is_valid)
            if (.not. is_valid) call write_parse_error_and_stop(subroutine_name, string)
            call create_valid_pattern(allow_int_min, allow_int_max, pattern, invalid_input_message)
        end subroutine check_args

        subroutine check_valid_input(char_idx)
            implicit none
            integer, intent(in) :: char_idx
            call check_substring(string(char_idx:char_idx), pattern, is_valid)
            if (.not. is_valid) then
                ! A number is NOT a integer or invalid input.
                if (rank == 0) print *, invalid_input_message, string(char_idx:)
                call write_parse_error_and_stop(subroutine_name, string)
            end if
        end subroutine check_valid_input

        subroutine read_rightnum(readnum)
            implicit none
            integer, intent(out) :: readnum
            read (string(rightnum_first_idx:), *, iostat=iostat) readnum
            if (iostat /= 0) then
                if (rank == 0) print *, "Can't get rightnum. string:", string, "rightnum", readnum
                call write_parse_error_and_stop(subroutine_name, string)
            end if
            call is_in_range_int(readnum, allow_int_min, allow_int_max, is_valid)
            if (.not. is_valid) then ! Rightnum is out of range
                write (min_str, *) allow_int_min; write (max_str, *) allow_int_max
                if (rank == 0) print *, "ERROR: rightnum is out of range,", &
                    "[", trim(adjustl(min_str)), ",", trim(adjustl(max_str)), "]"
                call write_parse_error_and_stop(subroutine_name, string)
            end if
        end subroutine read_rightnum

        subroutine search_rightnum_idx(right_idx)
            implicit none
            integer, intent(out) :: right_idx
            idx = verify(string(first_dot_index:), " ,.")
            if (idx == 0) call write_parse_error_and_stop(subroutine_name, string)

            right_idx = idx + first_dot_index - 1 ! Set the first index of the right num in string
            call check_valid_input(right_idx)
        end subroutine search_rightnum_idx

        subroutine read_leftnum(readnum)
            implicit none
            integer, intent(out) :: readnum
            read (string(leftnum_idx:first_dot_index - 1), *, iostat=iostat) readnum
            if (iostat /= 0) then
                if (rank == 0) print *, "Can't get leftnum. string:", string, "leftnum", readnum
                call write_parse_error_and_stop(subroutine_name, string)
            end if
            call is_in_range_int(readnum, allow_int_min, allow_int_max, is_valid)
            if (.not. is_valid) then ! The leftnum is out of range
                write (min_str, *) allow_int_min; write (max_str, *) allow_int_max
                if (rank == 0) print *, "ERROR: rightnum is out of range,", &
                    "[", trim(adjustl(min_str)), ",", trim(adjustl(max_str)), "]"
                call write_parse_error_and_stop(subroutine_name, string)
            end if
        end subroutine read_leftnum

        subroutine search_leftnum_idx(left_idx)
            integer, intent(out) :: left_idx
            left_idx = first_dot_index ! first dot '..' index in string
            do while (left_idx > 1)
                left_idx = left_idx - 1 ! (e.g. "2..5" -> "12..5")
                stat = verify(string(left_idx:first_dot_index), " ,") ! stat must be 1 or 2
                ! If stat is 2, we found the index of left num, so exit loop (e.g. string(left_idx:first_dot_index) = ",10.")
                if (stat == 2) then
                    left_idx = left_idx + 1  ! (e.g. string(left_idx:first_dot_index) = ",10." -> "10.")
                    exit
                else if (stat > 2 .or. stat <= 0) then
                    if (rank == 0) print *, "Can't get left num.string:", string, ",substring:", string(left_idx:first_dot_index)
                    call write_parse_error_and_stop(subroutine_name, string)
                end if
            end do
            call check_valid_input(left_idx)
        end subroutine search_leftnum_idx

        subroutine fill_numbers_to_list(left, right)
            implicit none
            integer, intent(in) :: left, right
            if (right < left) then
                ! right must be larger than or equal to left
                if (rank == 0) print *, "The specification of the range is invalid. left", left, "right", right
                call write_parse_error_and_stop(subroutine_name, string)
            end if
            ! Can fill numbers?
            if (size(list, 1) - filled_num < right - left + 1) then
                ! Can't fill numbers because the size of the list
                if (rank == 0) print *, "Can't fill range numbers because of the size of the list. size:", size(list, 1)
                call write_parse_error_and_stop(subroutine_name, string)
            end if

            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! Now we can fill the numbers
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            do idx = left, right
                filled_num = filled_num + 1
                list(filled_num) = idx
            end do
        end subroutine fill_numbers_to_list
    end subroutine parse_range_input_int

    subroutine read_subprograms(unit_num)
        use module_global_variables
        implicit none
        integer, intent(in) :: unit_num
        integer :: iostat
        character(len=max_str_length) :: input
        character(:), allocatable :: trim_input
        logical :: is_comment

        do while (.true.)
            read (unit_num, '(A)', iostat=iostat) input
            if (iostat /= 0) then
                if (rank == 0) print *, "ERROR: while reading subprograms, iostat = ", iostat, ", input =", input
                call stop_with_errorcode(iostat)
                call exit(iostat)
            end if

            call uppercase(input)
            call is_comment_line(input, is_comment)
            if (is_comment) cycle
            allocate (trim_input, source=trim(adjustl(input)))
            if (index(trim_input, ".") == 1) then
                ! If the input starts with a dot, it is the end of the subprograms.
                ! Need to reset the file pointer to the beginning of the line.
                backspace (unit_num)
                exit
            end if

            select case (trim_input)
            case ("CASCI")
                docasci = .true.
            case ("CASPT2")
                docaspt2 = .true.
            case ("IVO")
                doivo = .true.
            case default
                if (rank == 0) print *, "ERROR: Unknown subprogram: ", trim_input
                call stop_with_errorcode(1)
            end select
            deallocate (trim_input)
        end do

        if (doivo .and. (docasci .or. docaspt2)) then
            if (rank == 0) print *, "ERROR: ivo and casci or casci2 cannot be specified at the same time."
            call stop_with_errorcode(1)
        end if

    end subroutine read_subprograms

    subroutine read_caspt2_ciroots(unit_num)
        use module_global_variables
        use module_sort_swap, only: heapSort
        implicit none
        integer, intent(in) :: unit_num
        integer :: iostat, read_int, idx_filled, i, ciroots_idx, tmp_totsym, max_totsym
        logical :: is_comment
        integer :: tmp_ciroot_input(root_max), tmp_max_selectroots(totsym_max)
        integer, allocatable :: tmp_ciroots(:, :)
        character(len=max_str_length) :: input
        character(:), allocatable :: trim_input

        allocate (tmp_ciroots(ciroots_max, 2))
        tmp_ciroots(:, :) = 0
        tmp_max_selectroots(:) = 0
        ciroots_idx = 0
        max_totsym = 0
        do while (.true.)
            read (unit_num, '(A)', iostat=iostat) input
            if (iostat /= 0) then
                if (rank == 0) print *, "ERROR: while reading capst2_ciroots, iostat = ", iostat, ", input =", input
                call stop_with_errorcode(iostat)
                call exit(iostat)
            end if
            call is_comment_line(input, is_comment)
            if (is_comment) cycle

            allocate (trim_input, source=trim(adjustl(input)))
            if (index(trim_input, ".") == 1) then
                ! If the input starts with a dot, it is the end of the subprograms.
                ! Need to reset the file pointer to the beginning of the line.
                backspace (unit_num)
                exit
            end if
            call read_ciroots_line
            deallocate (trim_input)
            tmp_totsym = tmp_ciroot_input(1)
            max_totsym = max(max_totsym, tmp_totsym)
            call heapSort(list=tmp_ciroot_input(2:idx_filled), is_descending_order=.false.)
            tmp_max_selectroots(tmp_totsym) = max(tmp_ciroot_input(idx_filled), tmp_max_selectroots(tmp_totsym))

            ! Convert 1-dim list to 2-dim ciroots
            do i = 2, idx_filled
                ciroots_idx = ciroots_idx + 1
                tmp_ciroots(ciroots_idx, 1) = tmp_totsym ! total symmetry number
                tmp_ciroots(ciroots_idx, 2) = tmp_ciroot_input(i) ! selectroot
            end do

        end do

        allocate (caspt2_ciroots(ciroots_idx, 2))
        caspt2_ciroots = tmp_ciroots(1:ciroots_idx, :) ! Copy the tmp_ciroots to caspt2_ciroots
        deallocate (tmp_ciroots)

        allocate (max_selectroot_list(max_totsym))
        max_selectroot_list = tmp_max_selectroots(1:max_totsym) ! Copy the tmp_max_selectroots to max_selectroot_list
 
    contains
        subroutine read_ciroots_line
            implicit none
            integer :: idx_filled_by_range, idx_filled_int
            integer :: ciroots_range(root_max), ciroots_int(root_max)

            idx_filled_by_range = 0
            idx_filled_int = 0
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! First, we call this subroutine to detect the expression of the form a..b and expand it to a list.
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            call parse_range_input_int(string=trim_input, list=ciroots_range, filled_num=idx_filled_by_range, &
                                        allow_int_min=0, allow_int_max=root_max)

            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! Now, input string is seems to contain only space or , or integer
            ! So we parse integer from input string and expand it to a list
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            call parse_input_int(string=trim_input, list=ciroots_int, filled_num=idx_filled_int, &
                                allow_int_min=0, allow_int_max=root_max)

            idx_filled = idx_filled_by_range + idx_filled_int
            if (idx_filled - 1 > root_max) then ! first int in trim_input is totsym, therefore the number of root is idx_filled - 1
                if (rank == 0) then
                    print '(2(a,i0),2a)', "ERROR: .caspt2_ciroots max number of roots is ", root_max, &
                             " but selected ", idx_filled - 1, " roots. input:", input
                end if
                call stop_with_errorcode(1)                
            end if
            if (idx_filled < 2) then
                if (rank == 0) then
                    print *, "ERROR: .caspt2_ciroots must have at least 2 integers per line. input:", input
                    print *, "1st integer is the total symmetry number, 2nd integer and later are the selectroot."
                end if
                call stop_with_errorcode(1)
            end if
            tmp_ciroot_input(1:idx_filled_int) = ciroots_int(1:idx_filled_int)
            tmp_ciroot_input(idx_filled_int+1:idx_filled) = ciroots_range(1:idx_filled_by_range)
        end subroutine read_ciroots_line

    end subroutine read_caspt2_ciroots

    subroutine write_parse_error_and_stop(subroutine_name, input)
        implicit none
        character(len=*), intent(in) :: subroutine_name, input
        print *, "ERROR: Can't parse the input in ", subroutine_name, ", input:", input
        print *, "Stop the program. rank:", rank
        call stop_with_errorcode(1)
    end subroutine write_parse_error_and_stop

    subroutine check_range_allow_int(allow_int_min, allow_int_max, is_ok)
        implicit none
        integer, intent(in) :: allow_int_min, allow_int_max
        logical, intent(out) :: is_ok
        if (allow_int_max < allow_int_min) then
            if (rank == 0) print *, "ERROR: Allowed range of integer is invalid in parse_range_input_int.", &
                "MIN:", allow_int_min, "MAX:", allow_int_max
            is_ok = .false.
        else
            is_ok = .true.
        end if
    end subroutine check_range_allow_int

    subroutine check_substring(substring, string, is_substring_bool)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine returns whether varible "substring" is a substring of "string" or not
        ! (e.g.)  INPUT  : substring="ab" string="cdefgh"
        !         OUTPUT : is_substring_bool = .false.
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        implicit none
        character(*), intent(in) :: substring, string
        logical, intent(out) :: is_substring_bool
        if (index(string, substring) == 0) then
            is_substring_bool = .false.
        else
            is_substring_bool = .true.
        end if
    end subroutine check_substring

    subroutine is_in_range_int(num, num_min, num_max, is_in_range)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine returns whether varible "num" is in range of [num_min, num_max]
        ! (e.g.)  INPUT  : num = 10, num_min = -1, num_max = 11
        !         OUTPUT : is_in_range = .true.
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        implicit none
        integer, intent(in)     :: num, num_min, num_max
        logical, intent(out)    :: is_in_range
        if (num_min <= num .and. num <= num_max) then
            is_in_range = .true.
        else
            is_in_range = .false.
        end if
    end subroutine is_in_range_int

    subroutine create_valid_pattern(int_min, int_max, valid_pattern_string, invalid_message)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine returns a valid pattern string, valid_pattern_string, and a invalid_message in the event of an error.
        ! (e.g.)  INPUT  :  int_min = -11, int_max = -5
        !         OUTPUT :  valid_pattern_string = "-123456789"
        !                   invalid_message = "ERROR: Detected not minus numbers or Non-integer string."
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        implicit none
        integer, intent(in) :: int_min, int_max
        character(:), allocatable, intent(out) :: valid_pattern_string, invalid_message

        if (int_min > int_max) then
            if (rank == 0) print *, "ERROR: We want to create a allowed pattern of input, ", &
                "but the selected range of integer is invalid.", &
                "INT_MIN:", int_min, "INT_MAX", int_max
            if (rank == 0) print *, "Stop the program"
            call stop_with_errorcode(1)
        end if
        if (int_min < 0 .and. 0 <= int_max) then
            ! (e.g.) [-9, 10]
            valid_pattern_string = "-0123456789"
            invalid_message = "ERROR: Detected non integer string."
        elseif (int_min < 0) then
            ! int_min < 0 and int_max < 0
            ! (e.g.) [-10, -1]
            valid_pattern_string = "-123456789"
            invalid_message = "ERROR: Detected not minus numbers or non integer string."
        elseif (int_min == 0) then
            ! int_min = 0 and int_max >= 0
            ! (e.g.) [0, 100]
            valid_pattern_string = "0123456789"
            invalid_message = "ERROR: Detected minus numbers or non integer string."
        else
            ! int_min >= 0 and int_max >= 0
            ! (e.g.) [2, 7]
            valid_pattern_string = "123456789"
            invalid_message = "ERROR: Detected minus numbers or 0 or non integer string."
        end if
    end subroutine create_valid_pattern

    subroutine read_an_integer(unit_num, option_name, allowed_min_int, allowed_max_int, result_int)
        implicit none
        integer, intent(in) :: unit_num, allowed_min_int, allowed_max_int
        integer, intent(inout) :: result_int
        character(*), intent(in) :: option_name
        character(:), allocatable :: pattern, invalid_input_message
        logical :: is_comment, is_subst
        character(len=max_str_length) :: input
        call create_valid_pattern(allowed_min_int, allowed_max_int, pattern, invalid_input_message)
        do
            read (unit_num, '(a)') input
            call is_comment_line(input, is_comment)
            if (is_comment) cycle
            !  Is the input an integer and more than or equal to zero?
            call check_substring(input(1:1), pattern, is_subst)
            if (.not. is_subst) then
                if (rank == 0) then
                    print *, invalid_input_message, input
                    print *, 'invalidinput'
                end if
                call write_error_and_stop_read_an_integer
            end if
            read (input, *) result_int ! read an integer
            if (result_int < allowed_min_int .or. allowed_max_int < result_int) then
                if (rank == 0) then
                    print *, "ERROR: Your input is out of range. option_name: ", option_name
                    print *, "input:", result_int
                    print *, "minimum allowed value:", allowed_min_int
                    print *, "maximum allowed value:", allowed_max_int
                end if
                call write_error_and_stop_read_an_integer
            end if
            exit
        end do
    contains
        subroutine write_error_and_stop_read_an_integer
            implicit none
            if (rank == 0) then
                print *, "ERROR: Error in input, failed to read or validate an integer!!. Stop the program."
                print *, "input: ", trim(adjustl(input))
            end if
            call stop_with_errorcode(1)
        end subroutine write_error_and_stop_read_an_integer
    end subroutine read_an_integer

    subroutine read_a_string(unit_num, result_string)
        implicit none
        integer, intent(in) :: unit_num
        character(*), intent(inout) :: result_string
        logical :: is_comment
        character(len=max_str_length) :: input
        do
            read (unit_num, '(a)') input
            call is_comment_line(input, is_comment)
            if (is_comment) cycle
            read (input, *) result_string ! read a string
            exit
        end do
    end subroutine read_a_string

    subroutine is_comment_line(string, is_comment)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine returns trimmed input string and whether this input is comment line or not
        ! Comment line must include "!" or "#" in the first of character of input string except space
        ! (e.g.) INPUT  : string = "1,2,4..10,13!,17..20"
        !        OUTPUT : string = "1,2,4..10,13", is_comment = .false.
        !        INPUT2 : string = "#1,2,4..10,13"
        !        OUTPUT2: string = "#1,2,4..10,13", is_comment = .true.
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        implicit none
        character(*), intent(inout) :: string
        logical, intent(out) :: is_comment
        integer  :: comment_idx

        string = trim(adjustl(string)) ! (e.g.) "   2,3,4" => "2,3,4"

        ! Find the index of comment character (comment_idx = 0 if the comment character is not found)
        comment_idx = scan(string, '!#')
        if (verify(string, " ") == 0) then
            ! Empty line
            is_comment = .true.
        elseif (comment_idx == 1) then
            ! Comment line (e.g.) "!2,3,4"
            is_comment = .true.
        elseif (comment_idx == 0) then
            ! NOT Comment line (e.g.) "2,3,4"
            is_comment = .false.
        else
            ! NOT Comment line but include "!" or "#" (e.g.) "2,3,4 ! ras1"
            is_comment = .false.
            ! Trim before "!" or "#" (e.g.) "2,3,4 ! ras1" => "2,3,4 "
            string(:) = string(1:comment_idx - 1)
        end if

    end subroutine is_comment_line

    ! If the user enabled CASCI or CASPT2 subprograms, user must specify .caspt2_ciroots
    subroutine check_ciroots_set
        use module_global_variables
        implicit none
        if ((docasci .or. docaspt2) .and. .not. set_caspt2_ciroots) then
            if (rank == 0) then
                print *, "ERROR: you enabled CASCI or CASPT2 subprograms, but you didn't specify .caspt2_ciroots"
                print *, "Please specify .caspt2_ciroots in your input file."
            end if
            call stop_with_errorcode(1)
        end if
    end subroutine check_ciroots_set

    subroutine check_reqired_files_exist
        use module_global_variables
        implicit none
        logical :: is_exist, is_exist_mdcint

        ! Check the existence of MRCONEE file
        inquire (file="MRCONEE", exist=is_exist)
        if (.not. is_exist) then
            if (rank == 0) print *, "ERROR: MRCONEE file is required, but it is not found."
            call stop_with_errorcode(1)
        end if

        ! Check the existence of MDCINT or MDCINTNEW file
        inquire (file="MDCINT", exist=is_exist_mdcint)
        is_exist = is_exist_mdcint
        inquire (file="MDCINTNEW", exist=is_exist_mdcint)
        is_exist = or(is_exist, is_exist_mdcint)
        if (.not. is_exist) then
            if (rank == 0) print *, "ERROR: MDCINT or MDCINTNEW file is required, but it is not found."
            call stop_with_errorcode(1)
        end if

        if (doivo) then
            ! Check the existence of DFPCMO file
            inquire (file="DFPCMO", exist=is_exist)
            if (.not. is_exist) then
                if (rank == 0) print *, "ERROR: DFPCMO file is required, but it is not found."
                call stop_with_errorcode(1)
            end if
        end if

        if (.not. docasci .and. docaspt2) then
            ! Check the existence of CIDATA file
            inquire (file="CIDATA", exist=is_exist)
            if (.not. is_exist) then
                if (rank == 0) print *, "ERROR: CIDATA file is required, but it is not found."
                call stop_with_errorcode(1)
            end if
        end if
    end subroutine check_reqired_files_exist

    subroutine set_mdcint_scheme
        use module_global_variables, only: rank, is_scheme_set, mdcint_scheme, dirac_version, integrated_caspt2, &
                                           default_scheme_dirac22_or_earlier, default_scheme_dirac23_or_later
        implicit none

        ! If scheme option in active.inp is not set, set the default value.
        if (.not. is_scheme_set) then
            if (dirac_version > 22 .or. integrated_caspt2) then
                mdcint_scheme = default_scheme_dirac23_or_later
            else
                mdcint_scheme = default_scheme_dirac22_or_earlier
            end if
        end if
        is_scheme_set = .true.
        ! Validate mdcint_scheme value
        if (mdcint_scheme < 1) then
            ! Invalid mdcint_scheme input
            if (rank == 0) then
                print *, "Faild to validate the mdcint scheme value."
                print *, "mdcint_scheme must be larger than 1, but actual value:", mdcint_scheme
                print *, "Exit the program..."
            end if
            call stop_with_errorcode(1)
        end if

    end subroutine set_mdcint_scheme

    subroutine check_ras_is_valid
        use module_global_variables
        implicit none
        integer :: idx
        logical :: electron_filled(ninact + nact + nsec)

        ! min_hole_ras1 can't be larger than ras1_size
        if (min_hole_ras1 > ras1_size) then
            ! ERROR: The number of minimum hole of ras1 is larger than the number of ras1, It is unavailable.
            if (rank == 0) then
                print *, "ERROR: The number of minholeras1 is larger than the number of ras1."
                print *, "The number of ras1:", ras1_size
                print *, "The number of minholeras1:", min_hole_ras1
                print *, "The number of minholeras1 you specified is impossible."
                print *, "Exit the program."
                call stop_with_errorcode(1)
            end if
        end if

        ! min_elec_ras3 can't be larger than ras3_size
        if (min_elec_ras3 > ras3_size) then
            ! ERROR: The number of minimum electron of ras3 is larger than the number of ras3, It is unavailable.
            if (rank == 0) then
                print *, "ERROR: The number of minelecras3 is larger than the number of ras3."
                print *, "The number of ras3:", ras3_size
                print *, "The number of minelecras3:", min_elec_ras3
                print *, "The number of minelecras3 you specified is impossible."
                print *, "Exit the program."
                call stop_with_errorcode(1)
            end if
        end if

        ! min_elec_ras3 can't be larger than nelec
        if (min_elec_ras3 > nelec) then
            ! ERROR: The number of minimum electron of ras3 is larger than the number of active electron, It is unavailable.
            if (rank == 0) then
                print *, "ERROR: The number of minelecras3 is larger than the number of active electron."
                print *, "The number of active electron:", nelec
                print *, "The number of minelecras3:", min_elec_ras3
                print *, "The number of minelecras3 you specified is impossible."
                print *, "Exit the program."
                call stop_with_errorcode(1)
            end if
        end if      

        ! min_elec_ras3 can't be larger than ras3_max_elec
        if (min_elec_ras3 > ras3_max_elec) then
            ! ERROR: The number of minelecras3 is larger than the number of ras3_max_elec, It is unavailable.
            if (rank == 0) then
                print *, "ERROR: The number of minelecras3 is larger than ras3_max_elec."
                print *, "The number of ras3_max_elec:", ras3_max_elec
                print *, "The number of minelecras3:", min_elec_ras3
                print *, "The number of minelecras3 you specified is impossible."
                print *, "Exit the program."
                call stop_with_errorcode(1)
            end if
        end if      

        ! ras3_max_elec can't be larager than ras3_size
        if (ras3_max_elec > ras3_size) then
            ! ERROR: The number of max electron of ras3 is larger than the number of ras3, It is unavailable.
            if (rank == 0) then
                print *, "ERROR: The max number of allowed electron in ras3 is larger than the number of ras3."
                print *, "The number of ras3:", ras3_size
                print *, "The max number of allowed electron in ras3:", ras3_max_elec
                print *, "The max number of allowed electron in ras3 you specified is impossible."
                print *, "Exit the program."
                call stop_with_errorcode(1)
            end if
        end if

        ! Initialization
        electron_filled(:) = .false.

        ! Check duplication of electrons (ras1, ras2, ras3)
        if (ras1_size /= 0) then
            do idx = 1, ras1_size ! ras1_size is the size of the list.
                if (electron_filled(ras1_list(idx))) then
                    ! ERROR: The same number of the electron have been selected
                    if (rank == 0) print *, "ERROR: The number of selected more than once is", ras1_list(idx)
                    call write_error_and_stop_check_ras_is_valid ! Error in input. Stop the Program
                end if
                electron_filled(ras1_list(idx)) = .true. ! Fill ras1_list(idx)
            end do
        end if
        if (ras2_size /= 0) then
            do idx = 1, ras2_size ! ras2_size is the size of the list.
                if (electron_filled(ras2_list(idx))) then
                    ! ERROR: The same number of the electron have been selected
                    if (rank == 0) print *, "ERROR: The number of selected more than once is", ras2_list(idx)
                    call write_error_and_stop_check_ras_is_valid ! Error in input. Stop the Program
                end if
                electron_filled(ras2_list(idx)) = .true. ! Fill ras2_list(idx)
            end do
        end if
        if (ras3_size /= 0) then
            do idx = 1, ras3_size ! ras3_size is the size of the list.
                if (electron_filled(ras3_list(idx))) then
                    ! ERROR: The same number of the electron have been selected
                    if (rank == 0) print *, "ERROR: The number of selected more than once is", ras3_list(idx)
                    call write_error_and_stop_check_ras_is_valid ! Error in input. Stop the Program
                end if
                electron_filled(ras3_list(idx)) = .true. ! Fill ras3_list(idx)
            end do
        end if

        ! Is the number of RAS equal to the number of active?
        if (count(electron_filled) /= nact) then
            ! Error in input. Stop the Program
            if (rank == 0) then
                print *, "ERROR: Your input is invalid because the number of RAS is not equal to the number of active."
                print *, "active : ", nact
                print *, "RAS    : ", count(electron_filled)
                print *, "Stop the program."
            end if
            call stop_with_errorcode(1)
        end if
    contains
        subroutine write_error_and_stop_check_ras_is_valid
            implicit none
            if (rank == 0) then
                print *, "ERROR: Your input is invalid because the same number of the electron have been selected ", &
                    "in the RAS more than once!"
                print *, "YOUR INPUT"
                print *, "RAS1 : ", ras1_list
                print *, "RAS2 : ", ras2_list
                print *, "RAS3 : ", ras3_list
                print *, "Stop the program."
            end if
            call stop_with_errorcode(1)
        end subroutine write_error_and_stop_check_ras_is_valid
    end subroutine check_ras_is_valid

    subroutine lowercase(string)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine returns the lowercase string
        ! (e.g.) INPUT  : string = "tHiS"
        !        OUTPUT : string = "this"
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        implicit none
        character(*), intent(inout) :: string
        integer :: offset, idx, A_iachar, Z_iachar, chr_iachar
        offset = iachar('a') - iachar('A') ! offset number to convert uppercase to lowercase
        A_iachar = iachar('A') ! Ascii code of "A"
        Z_iachar = iachar('Z') ! Ascii code of "Z"
        do idx = 1, len(string)
            chr_iachar = iachar(string(idx:idx))
            if (A_iachar <= chr_iachar .and. chr_iachar <= Z_iachar) then
                ! A <= string(idx:idx) <= Z -> a <= string(idx:idx) <= z
                chr_iachar = chr_iachar + offset
                string(idx:idx) = achar(chr_iachar)
            end if
        end do
    end subroutine lowercase

    subroutine uppercase(string)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine returns the uppercase string
        ! (e.g.) INPUT  : string = "tHiS"
        !        OUTPUT : string = "THIS"
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        implicit none
        character(*), intent(inout) :: string
        integer :: offset, idx, a_iachar, z_iachar, chr_iachar
        offset = iachar('A') - iachar('a') ! offset number to convert uppercase to uppercase
        a_iachar = iachar('a') ! Ascii code of "a"
        z_iachar = iachar('z') ! Ascii code of "z"
        do idx = 1, len(string)
            chr_iachar = iachar(string(idx:idx))
            if (A_iachar <= chr_iachar .and. chr_iachar <= Z_iachar) then
                ! a <= string(idx:idx) <= z -> A <= string(idx:idx) <= Z
                chr_iachar = chr_iachar + offset
                string(idx:idx) = achar(chr_iachar)
            end if
        end do
    end subroutine uppercase

end module read_input_module
