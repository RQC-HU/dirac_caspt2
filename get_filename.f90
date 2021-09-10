subroutine get_mdcint_filename
    use four_caspt2_module, only: rank, mdcint_filename, mdcintnew, mdcint_debug, mdcint_int
    implicit none
    character(50)               :: file_basename, mdcint_basename, mdcintnum
    ! Rename the MDCINT to open according to the process number.
    file_basename = "MDCINXXXX"
    if (rank == 0) then
        mdcint_filename = "MDCINT"
        mdcintnew = "MDCINTNEW"
        mdcint_debug = "MDCINT_debug"
        mdcint_int = "MDCINT_int"
    else
        mdcint_basename = "MDCINXXXX"
        write (mdcintnum, "(I3)") rank
        mdcint_filename = TRIM(mdcint_baseName)//TRIM(ADJUSTL(mdcintnum))
        mdcintnew = "MDCINTNEW"//TRIM(ADJUSTL(mdcintnum))
        mdcint_debug = "MDCINT_debug"//TRIM(ADJUSTL(mdcintnum))
        mdcint_int = "MDCINT_int"//TRIM(ADJUSTL(mdcintnum))
    end if
    write(*,*)"get filename : ",trim(mdcint_filename)," ",trim(mdcintnew)," ",trim(mdcint_debug)," ",trim(mdcint_int)
end subroutine get_mdcint_filename

