!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Subroutine kr_copy_mpi ! 2 Electorn Integrals In Mdcint
program kr_copy_mpi

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    Use four_caspt2_module
    ! Use four_caspt2_module, only: indmor, nmo, realonly
    ! use omp_lib
    Implicit None
    include 'mpif.h'

    Character*50 :: Filename

    Character  :: datex*10, timex*8

    integer :: nkr = 96, nz
    integer :: i0, mdcint, inz, nnz
    integer :: ikr, jkr
    integer :: ii, jj, kk, ll
    integer :: iikr, jjkr, kkkr, llkr, iii, jjj, kkk, lll
    integer, allocatable :: indk(:), indl(:), kr(:)
    double precision, allocatable  :: rklr(:), rkli(:)
    ! Iwamuro modify
    real    :: cutoff
    integer :: nnkr, ikr8, jkr8, iiit, jjjt, kkkt, lllt
    integer, allocatable :: kkr(:), indk8(:), indl8(:)
    real*8, allocatable  :: rklr8(:), rkli8(:)

    ! integer         :: i, loop = 0, omp_max, tid
    Character*50    :: fileBaseName, mdcintBaseName, mdcintNew, mdcint_debug, mdcint_int, mdcintNum
    integer         :: ierr, nprocs, rank, procs = 8 !! TODO MPI PROCS を動的に?設定する
    real            :: time_start, time_end
    INTEGER         :: index
    nmo = 192
    ! omp_max = omp_get_max_threads()
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_rank(MPI_COMM_WORLD, rank, ierr)
    Allocate (kr(-nmo/2:nmo/2))
    kr = 0
    if (rank == 0) then ! Process limits for output
        ! open (8, file="debug", form="formatted", status="unknown")
        ! write (8, *) "debug open: rank=", rank
        write (*, *) "This is the serial execution area."
        ! do index = -nmo/2, nmo/2
        !     kr(index) = index
        ! end do
        open (10, file="MDCINT", form="unformatted", status="unknown")
        read (10) datex, timex, nkr, (kr(i0), kr(-1*i0), i0=1, nkr)
        open (11, file="firstline_mdcint", form="formatted", status="unknown")
        write (11, *) datex, timex, nkr, (kr(i0), kr(-1*i0), i0=1, nkr)
        ! write (11, *) size(kr)
        close (11)
        close (10)
    end if
    write (*, *) "end if: rank=", rank

    ! Allocate (indk(nmo**2))
    ! Allocate (indl(nmo**2))
    ! Allocate (rklr(nmo**2))
    ! Allocate (rkli(nmo**2))
    ! Allocate (rklr8(nmo**2))
    ! Allocate (rkli8(nmo**2))
    ! Allocate (indk8(nmo**2))
    ! Allocate (indl8(nmo**2))
    ! Allocate (kkr(-nmo/2:nmo/2))
    ! call MPI_Barrier(MPI_COMM_WORLD);
    ! write (*, *) "MPI_Barrier: rank=", rank
    ! call MPI_Barrier(MPI_COMM_WORLD);
    call MPI_Bcast(kr(-nmo/2), nmo + 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(datex, sizeof(datex), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(timex, sizeof(timex), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    write (*, *) "allocate successed. rank=", rank
    write (*, *) "if ierr == 0, kr broadcast successed. ierr=", ierr, "rank=", rank
    fileBaseName = "MDCINXXXX"
    if (rank == 0) then
        Filename = "MDCINT"
        mdcintNew = "MDCINTNEW"
        mdcint_debug = "MDCINT_debug"
        mdcint_int = "MDCINT_int"
    else
        mdcintBaseName = "MDCINXXXX"
        write (mdcintNum, "(I3)") rank
        Filename = TRIM(mdcintBaseName)//TRIM(ADJUSTL(mdcintNum))
        mdcintNew = "MDCINTNEW"//TRIM(ADJUSTL(mdcintNum))
        mdcint_debug = "MDCINT_debug"//TRIM(ADJUSTL(mdcintNum))
        mdcint_int = "MDCINT_int"//TRIM(ADJUSTL(mdcintNum))
    end if
    write (*, *) "set mdcint valiables. rank=", rank
    open (rank + 100, file=Filename, form='unformatted', status='unknown')
    open (rank + 300, file=mdcint_debug, form='formatted', status='unknown')

    read (rank + 100)
    write (rank + 300, *) datex, timex, nkr, (kr(i0), kr(-1*i0), i0=1, nkr)
    ! write (rank + 300, *) (kr(i0), kr(-1*i0), i0=1, nkr)
    close (rank + 100)
    close (rank + 300)
! 1000 close (rank + 100)
!     close (rank + 200)
!     close (rank + 300)
    ! close(30)
    ! deallocate (indk)
    ! deallocate (indl)
    ! deallocate (rklr)
    ! deallocate (rkli)
    ! deallocate (kkr)
    ! deallocate (indk8)
    ! deallocate (indl8)
    ! deallocate (rklr8)
    ! deallocate (rkli8)
    ! call MPI_Barrier(MPI_COMM_WORLD);
    call MPI_FINALIZE(ierr)
    write (20, *) "1000 closed "//trim(Filename)
    ! end do
    deallocate (kr)
    ! time_end = mpi_wtime()
    ! if (rank == 0) then
    ! close (8)
    ! open (10, file="output", form="formatted", status="unknown")
    ! write (10, *) "It took ", time_end - time_start, " seconds."
    ! close (10)
    ! end if
! end Subroutine kr_copy_mpi
end program kr_copy_mpi
