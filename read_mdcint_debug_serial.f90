program read_mdcint_debug
    implicit none
    integer                 ::  i = 7, mdcint = 11, debug = 12, ikr, jkr, nz, inz, nkr, i0, filesdebug = 13
    integer, allocatable    ::  indk(:), indl(:), kr(:)
    double precision, allocatable :: rklr(:), rkli(:)
    character*50    ::  Filename, mdcintBaseName, mdcint_debug, mdcintNum
    character*100   ::  errmsg
    integer         ::  nmo = 192

    inz = 1
    Allocate (kr(-nmo/2:nmo/2))
    Allocate (indk(nmo**2))
    Allocate (indl(nmo**2))
    Allocate (rklr(nmo**2))
    Allocate (rkli(nmo**2))
    mdcintBaseName = "MDCINXXXX"
    if (i == 1) then
        Filename = "MDCINT"
        ! mdcintNew = "MDCINTNEW"
        mdcint_debug = "MDCINT_debug_only"
        ! mdcint_int = "MDCINT_int"
    else
        write (mdcintNum, "(I3)") i - 1
        Filename = trim(mdcintBaseName)//trim(adjustl(mdcintNum))
        ! mdcintNew = "MDCINTNEW"
        mdcint_debug = "MDCINT_debug_only"//trim(adjustl(mdcintNum))
        ! mdcint_int = "MDCINT_int"
    end if
    open (mdcint, file=Filename, form="unformatted", status="unknown")
    open (debug, file=mdcint_debug, form="formatted", status="unknown")
    open (filesdebug, file="mdcintfiles_debug", form="formatted", status="unknown")
    write (filesdebug, *) Filename, mdcint_debug
    ! close(100)
    read (mdcint)
    ! read(mdcint, end=100, err=110, iomsg=errmsg) ikr,jkr, nz, (indk(inz),indl(inz), rklr(inz),rkli(inz), inz=1,nz)
    do
        read (mdcint, end=100, err=110, iomsg=errmsg) ikr, jkr, nz, &
                    (indk(inz), indl(inz), inz=1, nz), &
                    (rklr(inz), inz=1, nz)
        if (ikr == 0) then
            write (debug, "(3I20)") ikr, jkr, nz
            go to 100
        end if
        do inz = 1, nz
            write (debug, "(5I20,E32.16)") ikr, jkr, nz, indk(inz), indl(inz), rklr(inz)
        end do
    end do
100 write (filesdebug, *) "Read MDCINT"//trim(adjustl(mdcintNum))//" END"
    deallocate (kr)
    deallocate (indk)
    deallocate (indl)
    deallocate (rklr)
    deallocate (rkli)
    close (mdcint)
    close (debug)
    close (filesdebug)
110 write (filesdebug, *) "ERR : ", trim(errmsg)
    deallocate (kr)
    deallocate (indk)
    deallocate (indl)
    deallocate (rklr)
    deallocate (rkli)
    close (mdcint)
    close (debug)
    close (filesdebug)
end program read_mdcint_debug
