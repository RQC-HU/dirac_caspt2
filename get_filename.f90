subroutine get_mdcint_filename
    use four_caspt2_module, only: rank, mdcint_filename, mdcintnew, mdcint_debug, mdcint_int
    implicit none
    character(50)   :: mdcint_basename, mdcintnum, digit_x_padding
    ! Rename the MDCINT to open according to the process number.
    if (rank == 0) then
        mdcint_filename = "MDCINT"
        mdcintnew = "MDCINTNEW"
        mdcint_debug = "MDCINT_debug"
        mdcint_int = "MDCINT_int"
    else
        mdcint_basename = "MDCIN"
        if (rank >= 10000) then!! "ERROR": over five digit(can't assign)
            write (*, *) "ERROR: Can't assign MDCINT file to ranks of over five digits. rank:", rank
            stop
        else if (rank < 0) then !! "ERROR": minus number rank (can't assign)
            write (*, *) "ERROR: Can't assign MDCINT file to negative number of ranks. rank:", rank
            stop
        else if (rank < 10) then ! one digit (1~9)
            digit_x_padding = "XXXX"
        else if (rank < 100) then ! two digit (10~99)
            digit_x_padding = "XXX"
        else if (rank < 1000) then ! three digit (100~999)
            digit_x_padding = "XX"
        else if (rank < 10000) then ! four digit (1000~9999)
            digit_x_padding = "X"
        end if
        write (mdcintnum, "(I3)") rank
        mdcint_filename = TRIM(mdcint_baseName)//TRIM(ADJUSTL(digit_x_padding))//TRIM(ADJUSTL(mdcintnum))
        mdcintnew = "MDCINTNEW"//TRIM(ADJUSTL(mdcintnum))
        mdcint_debug = "MDCINT_debug"//TRIM(ADJUSTL(mdcintnum))
        mdcint_int = "MDCINT_int"//TRIM(ADJUSTL(mdcintnum))
    end if
    if (rank == 0) then
      write (3000, *) "get filename : ", trim(mdcint_filename), " ", trim(mdcintnew), " ", trim(mdcint_debug), " ", trim(mdcint_int)
    end if
end subroutine get_mdcint_filename
