module input_reader
    implicit none
    private
    public ras3_read
contains

    subroutine ras3_read
        use four_caspt2_module, only: rank, ras3_list, max_ras3_spinor_num
        implicit none
        integer :: idx, begin, digit
        integer, parameter :: length = 100
        character(length) :: line, ras3_chr, line_copy
        integer :: tmp_ras3(max_ras3_spinor_num)

        read (1, '(a)', err=10) line ! Read one line of active.inp
        line_copy = line
        idx = index(line, '..')
        print *, idx
        do while (idx /= 0)
            call parse_range_input_int(line, idx)
            idx = index(line, '..')
        end do

        idx = 1 ! The index of tmp_ras3
        ! The variable begin is the first index that does not contain space or , or ;

        begin = verify(line, ' ,') ! begin is 0 if all characters in line are space or , or ;
        do while (begin /= 0)
            print *, "begin", begin
            print *, 'line', line, line(begin:)
            read (line(begin:), *, err=10) tmp_ras3(idx) ! Read one of the ras3 value
            if (tmp_ras3(idx) <= 0) then ! Unexpected value error
                print *, "Unexpected value error!!, value:", tmp_ras3(idx)
                print *, "List  : ", tmp_ras3(1:idx)
                print *, "Input : ", line_copy
                print *, "Exit."
                stop
            end if
            write (ras3_chr, *) tmp_ras3(idx) ! ras3_chr is a string expression of tmp_ras3(idx)
            digit = len(trim(adjustl(ras3_chr))) ! The variable digit is the digit of tmp_ras3(idx) (e.g. -123->4, 10->2)
            line(begin:begin + digit) = "" ! Replace one of the ras3 value and separator to space (e.g. "10,3,5" -> "   3,5")
            print *, "After line", line
            begin = verify(line, ' ,-') ! Update the first index that does not contain space or , or ;
            idx = idx + 1 ! Count up the index of tmp_ras3
        end do
        idx = idx - 1
        allocate (ras3_list(idx)); Call memplus(KIND(ras3_list), SIZE(ras3_list), 1)
        ras3_list(:) = tmp_ras3(1:idx)
        print *, ras3_list
        goto 100
10      print *, "error in input, can't read ras3 value!!. Exit."
        stop
100     if (rank == 0) print *, "Read ras3 end"
    end subroutine ras3_read
    subroutine parse_range_input_int(line, first_dot_index)
        implicit none
        integer, intent(in) :: first_dot_index
        character(*), intent(inout) :: line
        character(30) :: right_str
        integer :: stat, rightnum_idx, leftnum_idx, leftnum, rightnum, rightnum_digit

        ! Get right num (e.g. "12..15" -> 15)
        rightnum_idx = verify(line(first_dot_index:), " ,.") ! Find the first index of right num
        rightnum_idx = rightnum_idx + first_dot_index - 1
        if (rightnum_idx == 0) goto 10 ! Right num is missing. Stop program
        read (line(rightnum_idx:), *, err=10) rightnum
        print *, "rightnum", rightnum
        write (right_str, *) rightnum
        rightnum_digit = len(trim(adjustl(right_str))) ! Get the digit of rightnum (e.g. -10 -> 3, 23 -> 2)

        ! Search the first index of the left num
        leftnum_idx = first_dot_index ! first dot '..' index in line
        do while (leftnum_idx > 0)
            leftnum_idx = leftnum_idx - 1 ! (e.g. "2." -> "12.")
            stat = verify(line(leftnum_idx:first_dot_index), " ,;") ! stat must be 1 or 2
            if (stat > 2 .or. stat <= 0) goto 10 ! Input Error. Stop program
            ! If stat is 2, exit loop (e.g. line(leftnum_idx:first_dot_index) = ",10.")
            if (stat == 2) then
                leftnum_idx = leftnum_idx + 1  ! (e.g. line(leftnum_idx:first_dot_index) = ",10." -> "10.")
                exit
            end if
        end do

        ! Now we can get the left num
        read (line(leftnum_idx:first_dot_index - 1), *) leftnum
        print *, "leftnum", leftnum
        line(leftnum_idx:rightnum_idx + rightnum_digit) = "" ! (e.g. "1, 2, 4-10, 13" -> "1, 2,      13" )
        print *, line
        ! Fill indices
        goto 100 ! End this subroutine
10      print *, "ERROR: Can't parse range input,", line, " Stop the program."; stop
100     continue
    end subroutine parse_range_input_int
end module input_reader
