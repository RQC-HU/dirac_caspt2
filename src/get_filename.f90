subroutine get_mdcint_filename(count)
    use four_caspt2_module, only: rank, nprocs, mdcint_filename, mdcintnew, mdcint_debug, mdcint_int
    implicit none
    character(50)   :: mdcint_basename, chr_rank, digit_x_padding
    integer, intent(in) :: count
    integer :: filename_idx
    ! Rename the MDCINT to open according to the process number.
    if (rank == 0 .and. count == 0) then
        mdcint_filename = "MDCINT"
        mdcintnew = "MDCINTNEW"
        mdcint_debug = "MDCINT_debug"
        mdcint_int = "MDCINT_int"
    else
        filename_idx = count*nprocs + rank
        mdcint_basename = "MDCIN"
        if (filename_idx >= 100000) then!! "ERROR": over six digit(can't assign)
            print *, "ERROR: Can't assign MDCINT file to ranks of over six digits. filename_idx:", filename_idx
            stop
        else if (filename_idx < 0) then !! "ERROR": minus number filename_idx (can't assign)
            print *, "ERROR: Can't assign MDCINT file to negative number of ranks. filename_idx:", filename_idx
            stop
        else if (filename_idx < 10) then ! one digit (1~9)
            digit_x_padding = "XXXX"
        else if (filename_idx < 100) then ! two digit (10~99)
            digit_x_padding = "XXX"
        else if (filename_idx < 1000) then ! three digit (100~999)
            digit_x_padding = "XX"
        else if (filename_idx < 10000) then ! four digit (1000~9999)
            digit_x_padding = "X"
        else if (filename_idx < 100000) then ! five digit (10000~99999)
            digit_x_padding = ""
        end if
        write (chr_rank, *) filename_idx
        mdcint_filename = TRIM(mdcint_baseName)//TRIM(ADJUSTL(digit_x_padding))//TRIM(ADJUSTL(chr_rank))
        if (count == 0) then
            mdcintnew = "MDCINTNEW"//TRIM(ADJUSTL(chr_rank))
            mdcint_debug = "MDCINT_debug"//TRIM(ADJUSTL(chr_rank))
            mdcint_int = "MDCINT_int"//TRIM(ADJUSTL(chr_rank))
        end if
    end if
    if (rank == 0) then
        print *, "get filename : ", trim(mdcint_filename), " ", &
            trim(mdcintnew), " ", trim(mdcint_debug), " ", trim(mdcint_int)
    end if
end subroutine get_mdcint_filename
subroutine get_subspace_filename
    use four_caspt2_module, only: rank, a1int, a2int, bint, c1int, c2int, c3int, &
                                  d1int, d2int, d3int, eint, fint, gint, hint
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
        write (chr_rank, *) rank
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
