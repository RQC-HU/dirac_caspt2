module module_error
    ! This module contains a subroutine to stop the program with a given error code.
    implicit none
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    private
    public :: stop_with_errorcode
contains
    subroutine stop_with_errorcode(errorcode)
        implicit none
        integer, intent(in) :: errorcode
        integer :: ierr
#ifdef HAVE_MPI
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
        print *, "Exit with error", errorcode
#ifdef INTEL
        call TRACEBACKQQ("Exit with error", errorcode)
#elif GNU
        call abort
#endif
        call exit(errorcode)
    end subroutine stop_with_errorcode
end module module_error
