program read_mdcint_debug
    implicit none
    include 'mpif.h'
    Character  :: datex*10, timex*8

    integer :: ikr, jkr, nz, inz, filesdebug = 1500, nkr, i0
    integer, allocatable :: indk(:), indl(:), kr(:)
    double precision, allocatable :: rklr(:), rkli(:)
    character*50        ::  Filename, mdcintBaseName, mdcint_debug, chr_mdcint
    integer             ::  nmo = 192
    real(8)             ::  time_start, time_end
    integer             ::  ierr, nprocs, rank, procs = 8 !! TODO MPI PROCS を動的に?設定する

    Allocate (kr(-nmo/2:nmo/2))
    kr = 0
    inz = 1
    ! add
    open (10, file="MDCINT", form="unformatted", status="unknown")
    read (10) datex, timex, nkr, (kr(i0), kr(-1*i0), i0=1, nkr)
    close (10)
    ! add end
    time_start = mpi_wtime()
    open (filesdebug, file="mdcintfiles_debug", form="formatted", status="unknown")
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_rank(MPI_COMM_WORLD, rank, ierr)
    Allocate (indk(nmo**2))
    Allocate (indl(nmo**2))
    Allocate (rklr(nmo**2))
    Allocate (rkli(nmo**2))
    mdcintBaseName = "MDCINXXXX"
    if (rank == 0) then
        Filename = "MDCINT"
        ! mdcintNew = "MDCINTNEW"
        mdcint_debug = "MDCINT_debug_only"
        ! mdcint_int = "MDCINT_int"
    else
        write (chr_mdcint, "(I3)") rank
        Filename = trim(mdcintBaseName)//trim(adjustl(chr_mdcint))
        ! mdcintNew = "MDCINTNEW"
        mdcint_debug = "MDCINT_debug_only"//trim(adjustl(chr_mdcint))
        ! mdcint_int = "MDCINT_int"
    end if
    write (200 + rank, *) 'debug', rank, trim(Filename), trim(mdcint_debug)

    open (rank + 10, file=Filename, form="unformatted", status="unknown")
    open (rank + 10 + nprocs, file=mdcint_debug, form="formatted", status="unknown")
    write (filesdebug, *) Filename, mdcint_debug
    read (rank + 10)
    do
        read (rank + 10, end=100) ikr, jkr, nz, &
            (indk(inz), indl(inz), inz=1, nz), &
            (rklr(inz), inz=1, nz)
        if (ikr == 0) then
            write (rank + 10 + nprocs, "(3I20)") ikr, jkr, nz
            write (filesdebug, "(3I20)") ikr, jkr, nz
            go to 100
        end if
        do inz = 1, nz
            write (rank + 10 + nprocs, "(5I20,E32.16)") ikr, jkr, nz, indk(inz), indl(inz), rklr(inz)
        end do
    end do
100 write (filesdebug, *) "Read MDCINT"//trim(adjustl(chr_mdcint))//" END"
    close (rank + 10)
    close (rank + 10 + nprocs)
    deallocate (indk)
    deallocate (indl)
    deallocate (rklr)
    deallocate (rkli)
    call MPI_FINALIZE(ierr)
    deallocate (kr)
    time_end = mpi_wtime()
    write (filesdebug, *) "It took ", time_end - time_start, " seconds.(RANK:", rank, ")"
    close (filesdebug)

end program read_mdcint_debug
