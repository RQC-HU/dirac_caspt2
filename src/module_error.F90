module module_error
    ! This module contains a subroutine to stop the program with a given error code.
    use, intrinsic :: iso_fortran_env, only: int32, int64
    implicit none
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    private
    public :: stop_with_errorcode
    interface stop_with_errorcode
        procedure stop_with_errorcode_i4
        procedure stop_with_errorcode_i8
    end interface stop_with_errorcode
contains
    subroutine stop_with_errorcode_i4(errorcode)
        implicit none
        integer(kind=int32), intent(in) :: errorcode
        integer(kind=int32) :: ierr
        ierr = 0
#ifdef INTEL
        call TRACEBACKQQ("Exit with error", errorcode)
#elif GNU
        call backtrace()
#endif
        call exit(errorcode)
    end subroutine stop_with_errorcode_i4

    subroutine stop_with_errorcode_i8(errorcode)
        implicit none
        integer(kind=int64), intent(in) :: errorcode
        integer(kind=int64) :: ierr
        ierr = 0
#ifdef INTEL
        call TRACEBACKQQ("Exit with error", errorcode)
#elif GNU
        call backtrace()
#endif
        call exit(int(errorcode, kind=int32))
    end subroutine stop_with_errorcode_i8
end module module_error
