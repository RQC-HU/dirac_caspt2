module input_reader
!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
! input_reader
! Copyright (c) by the authors of rel-caspt2.
! Author K.Noda
!
! This is a utility module that interpret and parse input strings.
!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!

    implicit none
    private
    public ras3_read, is_substring
    interface is_in_range_number
        module procedure is_in_range_int, is_in_range_real
    end interface is_in_range_number
contains

    subroutine ras3_read
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        ! This subroutine returns RAS3 list from the user input
        ! (e.g.) INPUT  : string = "1,2,4..10,13,17..20"
        !        OUTPUT : ras3_list = [1,2,4,5,6,7,8,9,10,13,17,18,19,20], (ras3_list is a global list in four_caspt2_module)
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        use four_caspt2_module, only: rank, ras3_list, max_ras3_spinor_num
        implicit none
        integer :: begin, digit
        integer, parameter :: max_str_length = 100
        character(max_str_length) :: string, ras3_chr, string_copy
        integer :: tmp_ras3(max_ras3_spinor_num), idx_filled

        read (1, '(a)', err=10) string ! Read a line of active.inp
        string_copy = string
        idx_filled = 0

        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        !  Parse the expressions of the form a..b in the input string and expands it to a list.
        !  (e.g.) INPUT  : string = "1,3,5..8,10", tmp_ras3 = [0,0,...,0],     idx_filled = 0
        !         OUTPUT : string = "1,3,    ,10", tmp_ras3 = [5,6,7,8,...,0], idx_filled = 4
        !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
        call parse_range_input_int(string, tmp_ras3, idx_filled, 0, 10**9)

        ! The variable begin is the first index that does not contain space or , or ;
        begin = verify(string, ' ,') ! begin is 0 if all characters in string are space or , or ;
        do while (begin /= 0)
            idx_filled = idx_filled + 1 ! Count up the index of tmp_ras3
            read (string(begin:), *, err=9) tmp_ras3(idx_filled) ! Read one of the ras3 value
            if (tmp_ras3(idx_filled) <= 0) then ! Unexpected value error
                print *, "ERROR: Unexpected value error!!, value:", tmp_ras3(idx_filled)
                print *, "List  : ", tmp_ras3(1:idx_filled)
                print *, "Input : ", string_copy
                print *, "Exit with an error."
                stop
            end if
            write (ras3_chr, *) tmp_ras3(idx_filled) ! ras3_chr is a string expression of tmp_ras3(idx_filled)
            digit = len(trim(adjustl(ras3_chr))) ! The variable digit is the digit of tmp_ras3(idx_filled) (e.g. -123->4, 10->2)
            string(begin:begin + digit - 1) = "" ! Replace one of the ras3 value and separator to space (e.g. "10,3,5" -> "   3,5")
            begin = verify(string, ' ,') ! Update the first index that does not contain space or , or ;
        end do
        allocate (ras3_list(idx_filled)); Call memplus(KIND(ras3_list), SIZE(ras3_list), 1)
        ras3_list(:) = tmp_ras3(1:idx_filled)
        print *, "ras3_list", ras3_list
        goto 100 ! Read the numbers properly
9       print *, "ERROR: Error in the section in reading the number, ", string(begin:)
10      print *, "ERROR: Error in input, can't read ras3 value!!. Stop the program."
        stop
100     if (rank == 0) print *, "Read ras3 end"
    end subroutine ras3_read

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
            print *, "LEFTNUM IDX:", leftnum_idx
            do while (leftnum_idx > 0)
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

end module input_reader
