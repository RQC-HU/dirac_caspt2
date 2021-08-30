program mpihello
    implicit none
    include 'mpif.h'
    integer         ::  nprocs, rank, ierr
    character(10)   ::  filename, chr_rank
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    write (chr_rank, "(I3)") rank
    filename = "output"//trim(adjustl(chr_rank))
    open (rank + 10, file=filename, form='formatted', status='unknown')
    write (rank + 10, '(a,2i8)') "Hello world Fortran", rank, nprocs
    close (rank + 10)
    call MPI_FINALIZE(ierr)
end program mpihello
! program main
! implicit none
! include 'mpif.h'

! integer i, size, rank, namelen, ierr
! character (len=MPI_MAX_PROCESSOR_NAME) :: name
! integer stat(MPI_STATUS_SIZE)

! call MPI_INIT (ierr)

! call MPI_COMM_SIZE (MPI_COMM_WORLD, size, ierr)
! call MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierr)
! call MPI_GET_PROCESSOR_NAME (name, namelen, ierr)

! if (rank.eq.0) then

!     print *, 'Hello world: rank ', rank, ' of ', size, ' running on ', name
!     open(10,file='output',form='formatted',status='unknown')
!     do i = 1, size - 1
!         call MPI_RECV (rank, 1, MPI_INTEGER, i, 1, MPI_COMM_WORLD, stat, ierr)
!         call MPI_RECV (size, 1, MPI_INTEGER, i, 1, MPI_COMM_WORLD, stat, ierr)
!         call MPI_RECV (namelen, 1, MPI_INTEGER, i, 1, MPI_COMM_WORLD, stat, ierr)
!         name = ''
!         call MPI_RECV (name, namelen, MPI_CHARACTER, i, 1, MPI_COMM_WORLD, stat, ierr)
!         write(10,*) 'Hello world: rank ', rank, ' of ', size, ' running on ', name
!         ! print *, 'Hello world: rank ', rank, ' of ', size, ' running on ', name
!     enddo

! else

!     call MPI_SEND (rank, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, ierr)
!     call MPI_SEND (size, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, ierr)
!     call MPI_SEND (namelen, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, ierr)
!     call MPI_SEND (name, namelen, MPI_CHARACTER, 0, 1, MPI_COMM_WORLD, ierr)

! endif

! call MPI_FINALIZE (ierr)

! end
