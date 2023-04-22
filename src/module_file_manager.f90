module module_file_manager
!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
! module_file_manager
! Copyright (c) by the authors of rel-caspt2.
! Author K.Noda
!
! This is a utility module that manages the file unit number.
!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
    use module_error, only: stop_with_errorcode
    use read_input_module, only: lowercase
    implicit none

contains

    subroutine check_iostat(iostat, file, end_of_file_reached)
        implicit none
        integer, intent(in) :: iostat
        character(len=*), intent(in) :: file
        logical, intent(out) :: end_of_file_reached
        if (iostat == 0) then
            end_of_file_reached = .false.
        else if (iostat < 0) then
            print *, "END OF FILE: ", file
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

    subroutine open_file(unit, form, file, status, action)
        implicit none
        character(len=*), intent(in) :: form, file, status, action
        integer, intent(inout) :: unit
        character(:), allocatable :: file_status
        integer :: iostat
        call search_unused_file_unit(unit)
        allocate (file_status, source=trim(status))
        call lowercase(file_status)
        if (file_status /= 'old' .and. file_status /= 'new' .and. file_status /= 'replace') then
            print *, 'ERROR: file_status must be old, new or replace. file_status = ', file_status
            print *, 'Exiting...'
            call stop_with_errorcode(1)
        end if
        open (unit, form=form, file=file, status=status, iostat=iostat, action=action)
        call check_file_open(file, iostat, unit)
    end subroutine open_file
    subroutine check_action_type(action, file)
        implicit none
        character(len=*), intent(in) :: action, file
        if (action /= 'read' .and. action /= 'write' .and. action /= 'readwrite') then
            print *, 'ERROR: action must be read, write or readwrite. action = ', action
            print *, 'FILE NAME: ', file
            print *, 'Exiting...'
            call stop_with_errorcode(1)
        end if
    end subroutine check_action_type

    subroutine open_unformatted_file(unit, file, status, optional_action)
        implicit none
        character(len=*), intent(in), optional :: optional_action
        character(len=*), intent(in) :: file, status
        integer, intent(inout) :: unit
        character(:), allocatable :: actual_action, trimmed_action, form

        if (present(optional_action)) then
            allocate (trimmed_action, source=trim(optional_action))
            call lowercase(trimmed_action)
            call check_action_type(action=trimmed_action, file=file)
            allocate (actual_action, source=trimmed_action)
        else
            allocate (actual_action, source='readwrite')
        end if
        allocate (form, source='unformatted')
        call open_file(unit=unit, form=form, file=file, status=status, action=actual_action)
    end subroutine open_unformatted_file

    subroutine open_formatted_file(unit, file, status, optional_action)
        implicit none
        character(len=*), intent(in), optional :: optional_action
        character(len=*), intent(in) :: file, status
        integer, intent(inout) :: unit
        character(:), allocatable :: form, actual_action, trimmed_action

        if (present(optional_action)) then
            allocate (trimmed_action, source=trim(optional_action))
            call lowercase(trimmed_action)
            call check_action_type(action=trimmed_action, file=file)
            allocate (actual_action, source=trimmed_action)
        else
            allocate (actual_action, source='readwrite')
        end if
        allocate (form, source='formatted')
        call open_file(unit=unit, form=form, file=file, status=status, action=actual_action)
    end subroutine open_formatted_file
end module module_file_manager
