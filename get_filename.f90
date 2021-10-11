subroutine get_mdcint_filename
    use four_caspt2_module, only: rank, mdcint_filename, mdcintnew, mdcint_debug, mdcint_int, normaloutput
    implicit none
    character(50)   :: mdcint_basename, chr_rank, digit_x_padding
    ! Rename the MDCINT to open according to the process number.
    ! if (rank == 0) then
    !     mdcint_filename = "MDCINT"
    !     mdcintnew = "MDCINTNEW1"
    !     mdcint_debug = "MDCINT_debug1"
    !     mdcint_int = "MDCINT_int1"
    ! else if (rank == 1) then
    !     mdcint_filename = "MDCINT"
    !     mdcintnew = "MDCINTNEW"
    !     mdcint_debug = "MDCINT_debug"
    !     mdcint_int = "MDCINT_int"
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
        write (chr_rank, "(I3)") rank
        mdcint_filename = TRIM(mdcint_baseName)//TRIM(ADJUSTL(digit_x_padding))//TRIM(ADJUSTL(chr_rank))
        mdcintnew = "MDCINTNEW"//TRIM(ADJUSTL(chr_rank))
        mdcint_debug = "MDCINT_debug"//TRIM(ADJUSTL(chr_rank))
        mdcint_int = "MDCINT_int"//TRIM(ADJUSTL(chr_rank))
    end if
    if (rank == 0) then ! Process limits for output
        write (normaloutput, *) "get filename : ", trim(mdcint_filename), " ", &
            trim(mdcintnew), " ", trim(mdcint_debug), " ", trim(mdcint_int)
    end if
end subroutine get_mdcint_filename
subroutine get_subspace_filename
    use four_caspt2_module, only: rank, a1int, a2int, bint, c1int, c2int, c3int, &
                                  d1int, d2int, d3int, eint, fint, gint, hint, normaloutput
    implicit none
    character(50) :: chr_rank
    if (rank == 0) then
        a1int = "A1int"
        a2int = "A2int"
        bint = "Bint"
        c1int = "C1int"
        c2int = "C2int"
        c3int = "C3int"
        d1int = "D1int"
        d2int = "D2int"
        d3int = "D3int"
        eint = "Eint"
        fint = "Fint"
        gint = "Gint"
        hint = "Hint"
    else
        write (chr_rank, "(I3)") rank
        a1int = "A1int"//TRIM(ADJUSTL(chr_rank))
        a2int = "A2int"//TRIM(ADJUSTL(chr_rank))
        bint = "Bint"//TRIM(ADJUSTL(chr_rank))
        c1int = "C1int"//TRIM(ADJUSTL(chr_rank))
        c2int = "C2int"//TRIM(ADJUSTL(chr_rank))
        c3int = "C3int"//TRIM(ADJUSTL(chr_rank))
        d1int = "D1int"//TRIM(ADJUSTL(chr_rank))
        d2int = "D2int"//TRIM(ADJUSTL(chr_rank))
        d3int = "D3int"//TRIM(ADJUSTL(chr_rank))
        eint = "Eint"//TRIM(ADJUSTL(chr_rank))
        fint = "Fint"//TRIM(ADJUSTL(chr_rank))
        gint = "Gint"//TRIM(ADJUSTL(chr_rank))
        hint = "Hint"//TRIM(ADJUSTL(chr_rank))
    end if

end subroutine get_subspace_filename
