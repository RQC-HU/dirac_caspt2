module module_file_manager
!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
! module_file_manager
! Copyright (c) by the authors of DIRAC-CASPT2.
! Author K.Noda
!
! This is a utility module that manages the file unit number.
!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
    use module_error, only: stop_with_errorcode
    use read_input_module, only: lowercase
    implicit none

    private
    public check_iostat, is_eof, open_unformatted_file, open_formatted_file

contains

    logical function is_eof(unit, file, is_formatted)
        ! Check whether a unit file is EOF or not.
        integer, intent(in) :: unit
        character(len=*), intent(in) :: file
        logical, intent(in) :: is_formatted
        integer :: iostat
        character(10) :: test_str_for_eof

        if (is_formatted) then
            read (unit, '(A)', iostat=iostat) test_str_for_eof
        else
            read (unit, iostat=iostat) test_str_for_eof
        end if
        call check_iostat(iostat, file, is_eof)
        backspace (unit) ! Since read test_str_for_eof, backspace to return to the begginning of the line
    end function

    subroutine check_iostat(iostat, file, end_of_file_reached)
        ! Check the value of iostat.
        ! If iostat < 0, then the end of file is reached.
        ! If iostat == 0, then the end of file is not reached.
        ! If iostat > 0, then an error occured.
        implicit none
        integer, intent(in) :: iostat
        character(len=*), intent(in) :: file
        logical, intent(out) :: end_of_file_reached
        if (iostat == 0) then
            end_of_file_reached = .false.
        else if (iostat < 0) then
            end_of_file_reached = .true.
        else
            print *, "ERROR: Error occured while reading a file. file: ", file, " iostat: ", iostat
            print *, "EXIT PROGRAM"
            call stop_with_errorcode(iostat)
        end if
    end subroutine check_iostat

    subroutine search_unused_file_unit(unit_file_num)
        implicit none
        integer, intent(inout) :: unit_file_num
        logical :: opened
        ! unit_file_num must be >= 21
        unit_file_num = 21
        ! Search for unused file unit
        do
            inquire (unit_file_num, opened=opened)
            if (.not. opened) exit ! unit_file_num is unused, so we can use it
            unit_file_num = unit_file_num + 1 ! Increment unit_file_num if the previous one is used
        end do
    end subroutine search_unused_file_unit

    subroutine check_file_open(file, iostat, unit)
        implicit none
        character(len=*), intent(in) :: file
        integer, intent(in) :: iostat, unit
        if (iostat .ne. 0) then
            print *, 'ERROR: Failed to open ', file, ': iostat = ', iostat, ' unit = ', unit
            print *, 'Exiting...'
            call stop_with_errorcode(iostat)
        end if
    end subroutine check_file_open

    subroutine open_file(unit, form, file, status, optional_action, optional_position)
        implicit none
        character(len=*), intent(in) :: form, file, status
        character(len=*), intent(in), optional :: optional_action, optional_position
        integer, intent(inout) :: unit
        character(:), allocatable :: file_status, action, position
        integer :: iostat
        call search_unused_file_unit(unit)
        allocate (file_status, source=parse_file_status(status))
        allocate (action, source=parse_optional_action(optional_action))
        if (present(optional_position)) then
            allocate (position, source=trim(adjustl(optional_position)))
            call lowercase(position)
            open (unit, form=form, file=file, status=file_status, iostat=iostat, action=action, position=position)
        else
            open (unit, form=form, file=file, status=file_status, iostat=iostat, action=action)
        end if
        call check_file_open(file, iostat, unit)
    end subroutine open_file

    function parse_file_status(status) result(parsed_status)
        implicit none
        character(len=*), intent(in) :: status
        character(:), allocatable :: parsed_status
        allocate (parsed_status, source=trim(adjustl(status)))
        call lowercase(parsed_status)
    end function parse_file_status

    function parse_optional_action(optional_action) result(parsed_action)
        implicit none
        character(len=*), intent(in), optional :: optional_action
        character(:), allocatable :: parsed_action
        if (present(optional_action)) then
            allocate (parsed_action, source=trim(adjustl(optional_action)))
            call lowercase(parsed_action)
        else
            allocate (parsed_action, source='readwrite')
        end if
    end function parse_optional_action

    subroutine open_unformatted_file(unit, file, status, optional_action, optional_position)
        implicit none
        character(len=*), intent(in), optional :: optional_action, optional_position
        character(len=*), intent(in) :: file, status
        integer, intent(inout) :: unit
        character(:), allocatable :: form

        allocate (form, source='unformatted')
        call open_file(unit, form, file, status, optional_action, optional_position)
    end subroutine open_unformatted_file

    subroutine open_formatted_file(unit, file, status, optional_action, optional_position)
        implicit none
        character(len=*), intent(in), optional :: optional_action, optional_position
        character(len=*), intent(in) :: file, status
        integer, intent(inout) :: unit
        character(:), allocatable :: form, actual_action, trimmed_action

        allocate (form, source='formatted')
        call open_file(unit, form, file, status, optional_action, optional_position)
    end subroutine open_formatted_file
end module module_file_manager
