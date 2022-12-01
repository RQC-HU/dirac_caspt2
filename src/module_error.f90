module module_error
    use four_caspt2_module, only: ierr
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
#ifdef HAVE_MPI
        call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
#endif
        call exit(errorcode)
    end subroutine stop_with_errorcode
end module module_error

