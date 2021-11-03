program tab_num_test
    implicit none
    integer::i, mdcint = 11, debug = 12,  inz
    integer(4) :: ikr, jkr, nz
    integer(4), allocatable :: indk(:), indl(:)
    double precision, allocatable :: rklr(:), rkli(:)
    character*50::Filename, mdcintBaseName, mdcint_debug, mdcintNum
    character*100:: errmsg

    integer :: nmo = 192
    inz = 1
    ! do i = 1, 8
        Allocate(indk(nmo**2))
        Allocate(indl(nmo**2))
        Allocate(rklr(nmo**2))
        Allocate(rkli(nmo**2))
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
        ! close(100)
        read (mdcint)
        ! read(mdcint, end=100, err=110, iomsg=errmsg) ikr,jkr, nz, (indk(inz),indl(inz), rklr(inz),rkli(inz), inz=1,nz)
        ! do inz = 1, nz
        do
            read(mdcint,end=100, err=110, iomsg=errmsg) ikr,jkr, nz, (indk(inz),indl(inz), rklr(inz), inz=1,nz)
            do inz = 1, nz
                write(debug, "(5I20,E32.16)") ikr, jkr, nz, indk(inz), indl(inz), rklr(inz)
            end do
            ! if ( mod(j,4) == 0 ) then
            !     write(debug, "()")
            ! end if
        end do
        100 write(100,*) "Read MDCINT"//trim(adjustl(mdcintNum))//" END"
            deallocate(indk)
            deallocate(indl)
            deallocate(rklr)
            deallocate(rkli)
            close(mdcint)
            close(debug)
            close(100)
    ! end do
110 write(100,*) "ERR : ", trim(errmsg)
    deallocate(indk)
    deallocate(indl)
    deallocate(rklr)
    deallocate(rkli)
    close(mdcint)
    close(debug)
    close(100)
end program tab_num_test
