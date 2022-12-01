module module_error
    implicit none
    include 'mpif.h'
    private
    public :: stop_with_errorcode
contains
    subroutine stop_with_errorcode(errorcode)
        implicit none
        integer, intent(in) :: errorcode
# ifdef  MPI
        call MPI_ABORT(MPI_COMM_WORLD, errorcode, MPI_IERROR)
# endif
        call exit(errorcode)
    end subroutine stop_with_errorcode
end module module_error

