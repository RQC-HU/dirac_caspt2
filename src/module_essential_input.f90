module module_essential_input
    use module_error, only: stop_with_errorcode
    use module_global_variables, only: rank
    implicit none
    private
    public :: add_essential_input, update_esesential_input, &
              check_all_essential_inputs_specified, &
              essential_inputs
    type essential_input
        character(:), allocatable :: name
        logical :: is_specified
    end type essential_input
    type(essential_input), allocatable :: essential_inputs(:)
contains

    subroutine add_essential_input(name)
        implicit none
        character(*), intent(in) :: name
        integer :: idx
        type(essential_input), allocatable :: tmp(:)

        if (.not. allocated(essential_inputs)) then
            ! Allocate the first element
            idx = 1
            allocate (essential_inputs(idx))
        else
            ! Reallocate the array with new_size = current_size + 1
            idx = size(essential_inputs, 1) + 1
            allocate (tmp(idx))
            tmp(1:idx - 1) = essential_inputs
            ! essential_inputs is reallocated with size = idx
            call move_alloc(tmp, essential_inputs)
        end if
        ! Add the name and default is_specified value
        essential_inputs(idx)%name = trim(adjustl(name))
        essential_inputs(idx)%is_specified = .false.
    end subroutine add_essential_input

    subroutine update_esesential_input(name, is_specified)
        implicit none
        character(*), intent(in) :: name
        logical, intent(in) :: is_specified
        integer :: idx
        character(:), allocatable :: trimmed_name

        ! Search the name in essential_inputs
        trimmed_name = trim(adjustl(name))
        do idx = 1, size(essential_inputs, 1)
            if (essential_inputs(idx)%name == trimmed_name) then
                essential_inputs(idx)%is_specified = is_specified
                return
            end if
        end do
        ! If the name is not found, it is an error.
        if (rank == 0) print *, "ERROR: Unknown input: ", trimmed_name
        call stop_with_errorcode(1)
    end subroutine update_esesential_input

    subroutine check_all_essential_inputs_specified()
        implicit none
        integer :: idx
        do idx = 1, size(essential_inputs, 1)
            if (.not. essential_inputs(idx)%is_specified) then
                if (rank == 0) print *, "ERROR: You must specify a variable '", trim(essential_inputs(idx)%name), "' before end."
                call stop_with_errorcode(1)
            end if
        end do
    end subroutine check_all_essential_inputs_specified

end module module_essential_input
