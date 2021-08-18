program tab_num_test
    implicit none
    integer::i, j, readval, mdcint = 11, debug = 12
    character*50::Filename, mdcintBaseName, mdcint_debug, mdcintNum
    character*100:: errmsg
    do i = 1, 8
        mdcintBaseName = "MDCINXXXX"
        if ( i == 1 ) then
            Filename = "MDCINT"
            ! mdcintNew = "MDCINTNEW"
            mdcint_debug = "MDCINT_debug"
            ! mdcint_int = "MDCINT_int"
        else
            write(mdcintNum,"(I3)") i-1
            Filename = trim(mdcintBaseName)//trim(adjustl(mdcintNum))
            ! mdcintNew = "MDCINTNEW"
            mdcint_debug = "MDCINT_debug"//trim(adjustl(mdcintNum))
            ! mdcint_int = "MDCINT_int"
        end if
        open(mdcint, file=Filename, form="unformatted", status="unknown")
        open(debug, file=mdcint_debug, form="formatted", status="unknown")
        open(100, file="mdcintfiles_debug", form="formatted", status="unknown")
        write(100,*) Filename, mdcint_debug
        close(100)
        j = 1
        do
            read(mdcint, end=100, err=110, iomsg=errmsg) readval
            write(debug, "(I4)", advance="no") readval
            if ( mod(j,4) == 0 ) then
                write(debug, "()")
            end if
            j = j + 1
        end do
    100 print *, "Read MDCINT END"
    fdsfksdjfds
    end do
110 print *, "ERR : ", trim(errmsg)
end program tab_num_test
