program mpihello
    implicit none
    include 'mpif.h'
    integer         ::  nprocs, rank, ierr, i = 1, nmo = 96, index
    integer, allocatable ::  j(:)
    character(10)   ::  filename, chr_rank
    allocate (j(-nmo/2:nmo/2))

    open (20, file="jlist", form="formatted", status="unknown")
    write (*, *) "This is the serial execution area."
    do index = -nmo/2, nmo/2
        j(index) = index
        write (*, *) index, j(index)
    end do
    open (10, file="listsize", form="formatted", status="unknown")
    write (10, *) size(j)
    close (10)
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    ! if (rank == 0) then
    !     open (20, file="jlist", form="formatted", status="unknown")
    !     write (*, *) "This is the serial execution area."
    !     do index = -nmo/2, nmo/2
    !         j(index) = index
    !         write (*, *) index, j(index)
    !     end do
    !     open (10, file="listsize", form="formatted", status="unknown")
    !     write (10, *) size(j)
    !     close (10)
    ! end if
    call MPI_Bcast(j(-nmo/2), nmo + 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    write (chr_rank, "(I3)") rank
    filename = "output"//trim(adjustl(chr_rank))
    open (rank + 10, file=filename, form='formatted', status='unknown')
    write (rank + 10, *) "Hello world Fortran", rank, nprocs
    do index = -nmo/2, nmo/2
        write (rank + 10, *) index, j(index)
    end do
    close (rank + 10)
    deallocate (j)
    call MPI_FINALIZE(ierr)
    write (20, *) 'rank = ', rank
    close (20)
end program mpihello
! program mpihello
!     implicit none
!     integer :: ierr
!     ! call MPI_Init(ierr)
!     write (*, *) "Hello"
! end program mpihello
