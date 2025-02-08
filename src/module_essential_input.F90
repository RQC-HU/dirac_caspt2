module module_essential_input
    use module_error, only: stop_with_errorcode
    use module_global_variables, only: rank
    implicit none
    private
    public :: essential_inputs_container

    type essential_input
        character(:), allocatable :: name
        logical :: is_specified
    end type essential_input

    ! Use a derived type to encapsulate the essential inputs array
    type essential_inputs_container
        type(essential_input), allocatable :: inputs(:)
    contains
        procedure :: add_essential_input
        procedure :: update_essential_input
        procedure :: check_all_essential_inputs_specified
        procedure :: get_essential_input_idx
        procedure :: essential_input_is_specified
    end type essential_inputs_container

contains

    subroutine add_essential_input(this, name)
        implicit none
        class(essential_inputs_container), intent(inout) :: this
        character(*), intent(in) :: name
        integer :: idx
        type(essential_input), allocatable :: tmp(:)

        if (.not. allocated(this%inputs)) then
            ! Allocate the first element
            idx = 1
            allocate (this%inputs(idx))
        else
            ! Reallocate the array with new_size = current_size + 1
            idx = size(this%inputs, 1) + 1
            allocate (tmp(idx))
            tmp(1:idx - 1) = this%inputs
            ! this%inputs is reallocated with size = idx
            call move_alloc(tmp, this%inputs)
        end if
        ! Add the name and default is_specified value
        this%inputs(idx)%name = trim(adjustl(name))
        this%inputs(idx)%is_specified = .false.
    end subroutine add_essential_input

    subroutine update_essential_input(this, name, is_specified)
        implicit none
        class(essential_inputs_container), intent(inout) :: this
        character(*), intent(in) :: name
        logical, intent(in) :: is_specified
        integer :: idx
        character(:), allocatable :: trimmed_name

        ! Search the name in essential_inputs
        trimmed_name = trim(adjustl(name))
        do idx = 1, size(this%inputs, 1)
            if (this%inputs(idx)%name == trimmed_name) then
                this%inputs(idx)%is_specified = is_specified
                return
            end if
        end do
        ! If the name is not found, it is an error.
        if (rank == 0) print *, "ERROR: Unknown input: ", trimmed_name
        call stop_with_errorcode(1)
    end subroutine update_essential_input

    subroutine check_all_essential_inputs_specified(this)
        implicit none
        class(essential_inputs_container), intent(in) :: this
        integer :: idx
        do idx = 1, size(this%inputs, 1)
            if (.not. this%inputs(idx)%is_specified) then
                if (rank == 0) print *, "ERROR: You must specify a variable '", trim(this%inputs(idx)%name), "' before end."
                call stop_with_errorcode(1)
            end if
        end do
    end subroutine check_all_essential_inputs_specified

    function get_essential_input_idx(this, name) result(idx)
        implicit none
        class(essential_inputs_container), intent(in) :: this
        character(*), intent(in) :: name
        integer :: i, idx
        character(:), allocatable :: trimmed_name

        trimmed_name = trim(adjustl(name))
        do i = 1, size(this%inputs, 1)
            if (this%inputs(i)%name == trimmed_name) then
                idx = i
                return ! Found, return idx
            end if
        end do
        idx = -1 ! Not found
    end function get_essential_input_idx

    function essential_input_is_specified(this, name) result(is_specified)
        implicit none
        class(essential_inputs_container), intent(in) :: this
        character(*), intent(in) :: name
        logical :: is_specified
        integer :: idx

        idx = get_essential_input_idx(this, name)
        if (idx == -1) then
            if (rank == 0) print *, "ERROR: Unknown input: ", trim(adjustl(name))
            call stop_with_errorcode(1)
        end if
        is_specified = this%inputs(idx)%is_specified
    end function essential_input_is_specified

end module module_essential_input
