module read_input_module
!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
! read_input_module
! Copyright (c) by the authors of rel-caspt2.
! Author K.Noda
!
! This is a utility module that interpret and parse input strings.
!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!

    implicit none
    private
    public read_input, is_substring, ras_read, lowercase
    logical is_end
    interface is_in_range_number
        module procedure is_in_range_int, is_in_range_real
    end interface is_in_range_number
contains
    subroutine read_input
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine is the entry point to read active.inp
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!

        implicit none
        integer :: idx
        character(100) :: string
        character(11), allocatable :: essential_variable_names(:)
        logical :: is_comment, is_config_sufficient, is_variable_filled(11) = &
                   (/.false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false./)
        is_end = .false.
        essential_variable_names = &
            (/"ninact    ", "nact      ", "nsec      ", "nroot     ", "nelec     ", &
            &  "selectroot", "totsym    ", "ncore     ", "nbas      ", "ptgrp     ", "diracver  "/)
        open (5, file="active.inp", form="formatted")
        do while (.not. is_end)
            read (5, "(a)", end=10) string
            call is_comment_line(string, is_comment)
            if (is_comment) cycle ! Read the next line
            call check_input_type(string, is_variable_filled)
        end do
        is_config_sufficient = .true.
        do idx = 1, size(is_variable_filled, 1)
            if (.not. is_variable_filled(idx)) then
                print *, "ERROR: You must specify a variable "//trim(essential_variable_names(idx))//" before end."
                is_config_sufficient = .false.
            end if
        end do
        if (.not. is_config_sufficient) goto 11 ! Error in input. Stop the Program
        ! if (size(ras1_list, 1) > 0 .or. size(ras2_list, 1) > 0 .or. size(ras3_list, 1) > 0) call check_ras_is_valid
        close (5)
        return ! END SUBROUTINE
10      print *, "YOU NEED TO ADD 'end' in active.inp"
        stop
11      print *, "ERROR: Error in input, valiables you specified is insufficient!!. Stop the program."
        stop
    end subroutine read_input

    subroutine check_input_type(string, is_filled)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine recognize the type of input that follows from the next line
        ! and calls the subroutine that we must call
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        use four_caspt2_module
        implicit none
        character(*), intent(inout) :: string
        logical, intent(inout) :: is_filled(:)
        call lowercase(string)
        select case (trim(string))

        case ("ninact")
            call read_an_integer(0, 10**9, ninact)
            is_filled(1) = .true.

        case ("nact")
            call read_an_integer(0, 10**9, nact)
            is_filled(2) = .true.

        case ("nsec")
            call read_an_integer(0, 10**9, nsec)
            is_filled(3) = .true.

        case ("nelec")
            call read_an_integer(0, 10**9, nelec)
            is_filled(4) = .true.

        case ("nroot")
            call read_an_integer(0, 10**9, nroot)
            is_filled(5) = .true.

        case ("selectroot")
            call read_an_integer(0, 10**9, selectroot)
            is_filled(6) = .true.

        case ("totsym")
            call read_an_integer(0, 10**9, totsym)
            is_filled(7) = .true.

        case ("ncore")
            call read_an_integer(0, 10**9, ncore)
            is_filled(8) = .true.

        case ("nbas")
            call read_an_integer(0, 10**9, nbas)
            is_filled(9) = .true.

        case ("eshift")
            read (5, *) eshift

        case ("ptgrp")
            read (5, '(A)') ptgrp
            is_filled(10) = .true.

        case ("diracver")
            call read_an_integer(0, 10**9, dirac_version)
            is_filled(11) = .true.

        case ("ras1")
            call ras_read(ras1_list, 1)

        case ("ras2")
            call ras_read(ras2_list, 2)

        case ("ras3")
            call ras_read(ras3_list, 3)

        case ("end")
            is_end = .true.

        end select
    end subroutine check_input_type
    subroutine ras_read(ras_list, ras_num)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine returns RAS[1,2,3] list from the user input
        ! (e.g.) INPUT  : string = "1,2,4..10,13,17..20"
        !        OUTPUT : ras[1,2,3]_list = [1,2,4,5,6,7,8,9,10,13,17,18,19,20], (ras[1,2,3]_list is a global list in four_caspt2_module)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        use four_caspt2_module, only: max_ras_spinor_num
        use module_sort_swap, only: heapSort
        implicit none
        integer, allocatable, intent(inout) :: ras_list(:)
        integer, intent(in) :: ras_num
        character(100) :: tmp_ras_chr
        character(:), allocatable :: ras_chr
        integer, parameter :: max_str_length = 100
        character(max_str_length) :: string
        integer :: tmp_ras(max_ras_spinor_num), idx_filled

        ! Get the ras_num and store this to ras_chr
        write (tmp_ras_chr, *) ras_num
        ras_chr = trim(adjustl(tmp_ras_chr))

        read (5, '(a)', err=10) string ! Read a line of active.inp
        idx_filled = 0

        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        !  Parse the expressions of the form a..b in the input string and expands it to a list.
        !  (e.g.) INPUT  : string = "1,3,5..8,10", tmp_ras = [0,0,...,0],                idx_filled = 0
        !         OUTPUT : string = " , ,    ,  ", tmp_ras = [5,6,7,8,1,3,10,0,0,...,0], idx_filled = 7
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        call parse_input_string_to_int_list(string, tmp_ras, idx_filled, 0, 10**9)

        ! Does the input string contain at least one varible?
        if (idx_filled <= 0) then
            print *, "string:", string
            goto 10 ! Input Error. Stop program
        end if
        allocate (ras_list(idx_filled)); Call memplus(KIND(ras_list), SIZE(ras_list), 1)
        ras_list(:) = tmp_ras(1:idx_filled)
        call heapSort(ras_list, .false.) ! Sort the ras_list in ascending order (lower to higher)
        print *, "ras"//ras_chr//"_list", ras_list
        goto 100 ! Read the numbers properly
10      print *, "ERROR: Error in input, can't read ras"//ras_chr//" value!!. Stop the program."
        stop
100     return ! END SUBROUTINE NORMALLY
    end subroutine ras_read

    subroutine parse_input_string_to_int_list(string, list, filled_num, allow_int_min, allow_int_max)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine returns a list of integers
        ! It finds expressions of integer or the form a..b in the input string and expands it to a list.
        ! (e.g.) INPUT  : string = "1,2,4..10,13,17..20", list = [-3,-4,0,0,...,0],                             filled_num = 2
        !        OUTPUT : string = " , ,     ,  ,      ", list = [1,2,4,5,6,7,8,9,10,13,17,18,19,20,0,0,...,0], filled_num = 14
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        implicit none
        integer, intent(in) :: allow_int_min, allow_int_max ! Allow allow_int_min <= x <= allow_int_max
        character(*), intent(inout) :: string ! Input string
        integer, intent(inout) :: filled_num ! The number of numbers already filled in list
        integer, intent(inout) :: list(:) ! A integer list

        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! First, we call this subroutine to detect the expression of the form a..b and expand it to a list.
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        call parse_range_input_int(string, list, filled_num, allow_int_min, allow_int_max)

        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! Now, input string is seems to contain only space or , or integer
        ! So we parse integer from input string and expand it to a list
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        call parse_input_int(string, list, filled_num, allow_int_min, allow_int_max)

        goto 100 ! NORMAL END
100     continue
    end subroutine parse_input_string_to_int_list

    subroutine parse_input_int(string, list, filled_num, allow_int_min, allow_int_max)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine returns a list of integers
        ! It finds expressions of integer in the input string and expands it to a list.
        ! This subroutine can't detect expressions of the form a..b,
        ! so You "must" call this subroutine after call "parse_range_input_int" subroutine.
        ! (e.g.) INPUT  : string = "1,2,4..10,13,17..20", list = [1,2,0,0,...,0],                             filled_num = 2
        !        OUTPUT : string = "1,2,     ,13,      ", list = [1,2,4,5,6,7,8,9,10,17,18,19,20,0,0,...,0], filled_num = 13
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        implicit none
        integer, intent(in) :: allow_int_min, allow_int_max ! Allow allow_int_min <= x <= allow_int_max
        character(*), intent(inout) :: string ! Input string
        integer, intent(inout) :: filled_num ! The number of numbers already filled in list
        integer, intent(inout) :: list(:) ! A integer list
        character(30) :: min_str, max_str, read_int_str
        character(:), allocatable  :: pattern, invalid_input_message
        integer :: read_int, read_int_digit, idx
        logical :: is_valid

        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! Checks for positive and negative integers and sets the first character patten
        ! and the error message that allowed within that range.
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        call create_valid_pattern(allow_int_min, allow_int_max, pattern, invalid_input_message)

        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! Read the number in the input string and store it to the list
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!

        ! The variable idx is the first index that does not contain space or , or ;
        idx = verify(string, ' ,') ! idx is 0 if all characters in string are space or , or ;
        do while (idx /= 0)
            ! Check whether the strint(idx:idx) is valid
            call is_substring(string(idx:idx), pattern, is_valid)
            if (.not. is_valid) then
                ! Right number is NOT a integer or invalid input.
                print *, invalid_input_message, string(idx:)
                goto 10 ! Input Error. Stop program
            end if

            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! Now we can get the num (e.g. "12,15" -> 12)
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            read (string(idx:), *, err=9) read_int ! Read one of the ras3 value
            ! Check whether the read_int is in range [allow_int_min, allow_int_max]
            call is_in_range_number(read_int, allow_int_min, allow_int_max, is_valid)
            if (.not. is_valid) then
                ! The read_int is out of range [allow_int_min, allow_int_max]
                write (min_str, *) allow_int_min
                write (max_str, *) allow_int_max
                print *, "ERROR: read_int is out of range,", &
                    "[", trim(adjustl(min_str)), ",", trim(adjustl(max_str)), "]"
                goto 10 ! Input Error. Stop program
            end if

            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! Try to get the digit of read_int to rewrite the section of read_int as blank
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            write (read_int_str, *) read_int ! read_int_str is a string expression of read_int
            read_int_digit = len(trim(adjustl(read_int_str))) ! The variable read_int_digit is the read_int_digit of read_int (e.g. -123->4, 10->2)

            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! Check whether we can fill the numbers
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            filled_num = filled_num + 1 ! Count up the index of the list
            if (size(list, 1) < filled_num) then
                ! Can't fill numbers because the size of the list
                print *, "Can't fill range numbers because of the size of the list. size:", size(list, 1)
                goto 10 ! Input Error. Stop program
            end if

            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! Now we can store the number
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            list(filled_num) = read_int

            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! Rewrite the section we read as blank. (e.g.) "10,3,5" -> "  ,3,5"
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            string(idx:idx + read_int_digit - 1) = ""

            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! Find the next index that does not contain space or , or ;
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            idx = verify(string, ' ,')

        end do

        goto 100 ! NORMAL END
9       print *, "ERROR: Error in the section in reading the number, ", string(idx:); goto 10
10      print *, "ERROR: Can't parse the input in parse_input_int, input:", string, " Stop the program."
        stop
100     continue
    end subroutine parse_input_int

    subroutine parse_range_input_int(string, list, filled_num, allow_int_min, allow_int_max)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine returns a list of integers
        ! It finds expressions of the form a..b in the input string and expands it to a list.
        ! (e.g.) INPUT  : string = "1,2,4..10,13,17..20", list = [1,2,0,0,...,0],                             filled_num = 2
        !        OUTPUT : string = "1,2,     ,13,      ", list = [1,2,4,5,6,7,8,9,10,17,18,19,20,0,0,...,0], filled_num = 13
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        implicit none
        integer, intent(in) :: allow_int_min, allow_int_max ! Allow allow_int_min <= x <= allow_int_max
        character(*), intent(inout) :: string ! Input string
        integer, intent(inout) :: filled_num ! The number of numbers already filled in list
        integer, intent(inout) :: list(:) ! A integer list
        character(30) :: right_str, min_str, max_str
        character(:), allocatable  :: pattern, invalid_input_message
        integer :: first_dot_index, stat, rightnum_idx, leftnum_idx, leftnum, rightnum, rightnum_digit, idx
        logical :: is_valid

        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=
        ! Checks whether the range of allow_int is valid
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=
        if (allow_int_max < allow_int_min) then
            print *, "ERROR: Allowed range of integer is invalid in parse_range_input_int.", &
                "MIN:", allow_int_min, "MAX:", allow_int_max
            goto 10 ! Input Error. Stop program
        end if

        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! Checks for positive and negative integers and sets the first character patten
        ! and the error message that allowed within that range.
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        call create_valid_pattern(allow_int_min, allow_int_max, pattern, invalid_input_message)

        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! Find the first index of double dots ".."
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        first_dot_index = index(string, '..')

        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! Read the numbers of the left and right and fill in the list to the range of the readings
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        do while (first_dot_index /= 0) ! Find a ".." expression in the string

            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! Find the first index of the right num
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            rightnum_idx = verify(string(first_dot_index:), " ,.") ! Find the first index of the right num in string(first_dot_index:)
            if (rightnum_idx == 0) goto 10 ! Right num is missing. Stop program
            rightnum_idx = rightnum_idx + first_dot_index - 1 ! Set the first index of the right num in string
            ! Check whether the first character of the right num is valid
            call is_substring(string(rightnum_idx:rightnum_idx), pattern, is_valid)
            if (.not. is_valid) then
                ! Right number is NOT a integer or invalid input.
                print *, invalid_input_message, string(rightnum_idx:)
                goto 10 ! Input Error. Stop program
            end if

            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! Now we can get the right num (e.g. "12..15" -> 15)
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            read (string(rightnum_idx:), *, err=8) rightnum ! err=8 : Invalid input
            ! Check whether the rightnum is in range [allow_int_min, allow_int_max]
            call is_in_range_number(rightnum, allow_int_min, allow_int_max, is_valid)
            if (.not. is_valid) then
                ! The rightnum is out of range [allow_int_min, allow_int_max]
                write (min_str, *) allow_int_min
                write (max_str, *) allow_int_max
                print *, "ERROR: rightnum is out of range,", &
                    "[", trim(adjustl(min_str)), ",", trim(adjustl(max_str)), "]"
                goto 10 ! Input Error. Stop program
            end if
            write (right_str, *) rightnum
            rightnum_digit = len(trim(adjustl(right_str))) ! Get the digit of rightnum (e.g. -10 -> 3, 23 -> 2)

            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! Search the first index of the left num
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            leftnum_idx = first_dot_index ! first dot '..' index in string
            do while (leftnum_idx > 1)
                leftnum_idx = leftnum_idx - 1 ! (e.g. "2..5" -> "12..5")
                stat = verify(string(leftnum_idx:first_dot_index), " ,;") ! stat must be 1 or 2
                if (stat > 2 .or. stat <= 0) then
                    print *, "Can't get left num. substring:", string(leftnum_idx:first_dot_index)
                    goto 10 ! Input Error. Stop program
                end if
                ! If stat is 2, we found the index of left num, so exit loop (e.g. string(leftnum_idx:first_dot_index) = ",10.")
                if (stat == 2) then
                    leftnum_idx = leftnum_idx + 1  ! (e.g. string(leftnum_idx:first_dot_index) = ",10." -> "10.")
                    exit
                end if
            end do
            ! Check whether the first character of the left num is valid
            call is_substring(string(leftnum_idx:leftnum_idx), pattern, is_valid)
            if (.not. is_valid) then
                ! Right number is NOT a integer or invalid input.
                print *, invalid_input_message, string(leftnum_idx:)
                goto 10 ! Input Error. Stop program
            end if

            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! Now we can get the left num (e.g. "12..15" -> 12)
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            read (string(leftnum_idx:first_dot_index - 1), *, err=9) leftnum ! err=9 : Invalid input
            ! Check whether the rightnum is in range [allow_int_min, allow_int_max]
            call is_in_range_number(rightnum, allow_int_min, allow_int_max, is_valid)
            if (.not. is_valid) then
                ! The rightnum is out of range [allow_int_min, allow_int_max]
                write (min_str, *) allow_int_min
                write (max_str, *) allow_int_max
                print *, "ERROR: rightnum is out of range,", &
                    "[", trim(adjustl(min_str)), ",", trim(adjustl(max_str)), "]"
                goto 10 ! Input Error. Stop program
            end if

            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! Rewrite the section we read as blank. (e.g. "1, 2, 4..10, 13" -> "1, 2,      13" )
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            string(leftnum_idx:rightnum_idx + rightnum_digit - 1) = ""

            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! Check whether we can fill the numbers
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            if (rightnum < leftnum) then
                ! rightnum must be larger than or equal to leftnum
                print *, "The specification of the range is invalid. left", leftnum, "right", rightnum
                goto 10  ! Input Error. Stop program
            end if
            ! Can fill numbers?
            if (size(list, 1) - filled_num < rightnum - leftnum + 1) then
                ! Can't fill numbers because the size of the list
                print *, "Can't fill range numbers because of the size of the list. size:", size(list, 1)
                goto 10 ! Input Error. Stop program
            end if

            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! Now we can fill the numbers
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            do idx = leftnum, rightnum
                filled_num = filled_num + 1
                list(filled_num) = idx
            end do

            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            ! Find the next index of double dots ".."
            !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
            first_dot_index = index(string, '..')
        end do

        goto 100 ! End this subroutine
8       print *, "Can't get rightnum. string:", string, "rightnum", rightnum; goto 10 ! Stop program (error)
9       print *, "Can't get leftnum. string:", string, "leftnum:", leftnum; goto 10 ! Stop program (error)
10      print *, "ERROR: Can't parse the input in parse_range_input_int, input:", string, " Stop the program."
        stop ! Stop program (error)
100     continue ! Read the numbers properly
    end subroutine parse_range_input_int

    subroutine is_substring(substring, string, is_substring_bool)
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
    end subroutine is_substring

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

    subroutine is_in_range_real(num, num_min, num_max, is_in_range)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine returns whether varible "num" is in range of [num_min, num_max]
        ! (e.g.)  INPUT  : num = 10, num_min = -1, num_max = 11
        !         OUTPUT : is_in_range = .true.
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        implicit none
        real(8), intent(in)     :: num, num_min, num_max
        logical, intent(out)    :: is_in_range
        if (num_min <= num .and. num <= num_max) then
            is_in_range = .true.
        else
            is_in_range = .false.
        end if
    end subroutine is_in_range_real

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
            print *, "ERROR: We want to create a allowed pattern of input, but the selected range of integer is invalid.", &
                "INT_MIN:", int_min, "INT_MAX", int_max
            print *, "Stop the program"
            stop
        end if
        if (int_min < 0 .and. 0 <= int_max) then
            ! (e.g.) [-9, 10]
            valid_pattern_string = "-0123456789"
            invalid_message = "ERROR: Detected Non-integer string."
        elseif (int_min < 0) then
            ! int_min < 0 and int_max < 0
            ! (e.g.) [-10, -1]
            valid_pattern_string = "-123456789"
            invalid_message = "ERROR: Detected not minus numbers or Non-integer string."
        elseif (int_min == 0) then
            ! int_min = 0 and int_max >= 0
            ! (e.g.) [0, 100]
            valid_pattern_string = "0123456789"
            invalid_message = "ERROR: Detected minus numbers or Non-Integer string."
        else
            ! int_min >= 0 and int_max >= 0
            ! (e.g.) [2, 7]
            valid_pattern_string = "123456789"
            invalid_message = "ERROR: Detected minus numbers or 0 or Non-Integer string."
        end if
    end subroutine create_valid_pattern
    subroutine read_an_integer(allowed_min_int, allowed_max_int, result_int)
        implicit none
        integer, intent(in) :: allowed_min_int, allowed_max_int
        integer, intent(inout) :: result_int
        character(:), allocatable :: pattern, invalid_input_message
        logical :: is_comment, is_subst
        character(100) :: input
        call create_valid_pattern(allowed_min_int, allowed_max_int, pattern, invalid_input_message)
        do
            read (5, '(a)') input
            call is_comment_line(input, is_comment)
            if (is_comment) cycle ! Go to the next line
            !  Is the input an integer and more than or equal to zero?
            call is_substring(input(1:1), pattern, is_subst)
            if (.not. is_subst) then
                print *, invalid_input_message, input
                print *, 'invalidinput'
                goto 10
            end if
            read (input, '(I)') result_int ! read an integer
            exit ! EXIT LOOP
        end do
        return ! END SUBROUTINE
10      print *, "ERROR: Error in input, can't read a integer value!!. Stop the program."
        print *, "input: ", input
        stop
    end subroutine read_an_integer

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

        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! Trim string
        ! (e.g.) "   2,3,4" => "2,3,4"
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        string = trim(string)
        comment_idx = scan(string, '!#')
        if (comment_idx == 1) then
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

    subroutine check_ras_is_valid
        use four_caspt2_module, only: ras1_list, ras2_list, ras3_list, ninact, nact, nsec
        implicit none
        integer :: idx
        logical :: electron_filled(ninact + nact + nsec)
        electron_filled(:) = .false.
        do idx = 1, size(ras1_list, 1) ! size(ras1_list,1) is the size of the list.
            if (electron_filled(ras1_list(idx))) then
                ! ERROR: The same number of the electron have been selected
                print *, "ERROR: The number of selected more than once is", ras1_list(idx)
                goto 10 ! Error in input. Stop the Program
            end if
            electron_filled(ras1_list(idx)) = .true. ! Fill ras1_list(idx)
        end do
        do idx = 1, size(ras2_list, 1) ! size(ras2_list,1) is the size of the list.
            if (electron_filled(ras2_list(idx))) then
                ! ERROR: The same number of the electron have been selected
                print *, "ERROR: The number of selected more than once is", ras2_list(idx)
                goto 10 ! Error in input. Stop the Program
            end if
            electron_filled(ras2_list(idx)) = .true. ! Fill ras2_list(idx)
        end do
        do idx = 1, size(ras3_list, 1) ! size(ras3_list,1) is the size of the list.
            if (electron_filled(ras3_list(idx))) then
                ! ERROR: The same number of the electron have been selected
                print *, "ERROR: The number of selected more than once is", ras3_list(idx)
                goto 10 ! Error in input. Stop the Program
            end if
            electron_filled(ras3_list(idx)) = .true. ! Fill ras3_list(idx)
        end do

        ! Is the number of RAS equal to the number of active?
        if (count(electron_filled) /= nact) then
            goto 11 ! Error in input. Stop the Program
        end if

        return ! END NORMALLY
10      print *, "ERROR: Your input is invalid because the same number of the electron have been selected in the RAS" &
            //" more than once!"
        print *, "YOUR INPUT"
        print *, "RAS1 : ", ras1_list
        print *, "RAS2 : ", ras2_list
        print *, "RAS3 : ", ras3_list
        print *, "Stop the program."
        stop
11      print *, "ERROR: Your input is invalid because the number of RAS is not equal to the number of active."
        print *, "active : ", nact
        print *, "RAS    : ", count(electron_filled)
        print *, "Stop the program."
        stop
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
end module read_input_module