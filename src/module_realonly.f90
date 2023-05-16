module module_realonly
    implicit none
    private
    public :: realonly_type, check_realonly

    type :: realonly_type
        logical, private :: is_real
    contains
        private
        procedure, public :: set_is_realonly
        procedure, public :: is_realonly
    end type realonly_type

    type(realonly_type), public :: realonly
contains

    subroutine set_is_realonly(this, is_real)
        class(realonly_type), intent(inout) :: this
        logical, intent(in) :: is_real
        this%is_real = is_real
    end subroutine set_is_realonly

    function is_realonly(this) result(is_real)
        class(realonly_type), intent(in) :: this
        logical :: is_real
        is_real = this%is_real
    end function is_realonly

    subroutine check_realonly()
        use module_global_variables
        use module_error
        use module_file_manager

        integer :: unit_mdcint
        integer :: i, j, nz, inz, iostat
        integer, allocatable :: k(:), l(:)
        real(8), allocatable :: rklr(:), rkli(:)
        character(6), parameter :: filename = "MDCINT"

        allocate (k(nmo**2), l(nmo**2), rklr(nmo**2), rkli(nmo**2))

        call open_unformatted_file(unit_mdcint, filename, "old")
        rewind (unit_mdcint)
        read (unit_mdcint) ! Skip header
        read (unit_mdcint, iostat=iostat) i, j, nz, (k(inz), l(inz), inz=1, nz), (rklr(inz), rkli(inz), inz=1, nz)
        if (iostat == 0) then ! Complex
            call set_is_realonly(realonly, .false.)
        else ! Realonly or Error
            rewind (unit_mdcint) ! Go back to the beginning of the file
            read (unit_mdcint) ! Skip header
            read (unit_mdcint, iostat=iostat) i, j, nz, (k(inz), l(inz), inz=1, nz), (rklr(inz), inz=1, nz)
            if (iostat == 0) then ! Realonly
                call set_is_realonly(realonly, .true.)
            else ! Error
                if (rank == 0) print *, "Format of MDCINT is different from the one expected by this program."
                call stop_with_errorcode(iostat)
            end if
        end if
        if (rank == 0) print *, "MDCINT realonly = ", realonly%is_realonly()
        deallocate (k, l, rklr, rkli)
        close (unit_mdcint)
    end subroutine check_realonly

end module module_realonly
