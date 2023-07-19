module module_error
    ! This module contains a subroutine to stop the program with a given error code.
    implicit none
    private
    public :: stop_with_errorcode
contains
    subroutine stop_with_errorcode(errorcode)
        implicit none
        integer, intent(in) :: errorcode
        print *, "Exit with error", errorcode
#ifdef INTEL
        call TRACEBACKQQ("Exit with error", errorcode)
#elif GNU
        call abort
#endif
        call exit(errorcode)
    end subroutine stop_with_errorcode
end module module_error
