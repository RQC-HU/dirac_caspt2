module module_file_manager
!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
! module_file_manager
! Copyright (c) by the authors of rel-caspt2.
! Author K.Noda
!
! This is a utility module that manages the file unit number.
!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
    use read_input_module, only: lowercase
    implicit none

contains
    subroutine search_unused_file_unit(file_unit_number)
        implicit none
        integer, intent(inout) :: file_unit_number
        logical :: opened
        ! file_unit_number must be >= 20
        if (file_unit_number < 20) then
            file_unit_number = 20
        end if
        ! Search for unused file unit
        do
            file_unit_number = file_unit_number + 1
            inquire (file_unit_number, opened=opened)
            if (.not. opened) exit ! file_unit_number is unused, so we can use it
        end do
    end subroutine search_unused_file_unit

    subroutine check_file_open(file, iostat)
        implicit none
        character(len=*), intent(in) :: file
        integer, intent(in) :: iostat
        if (iostat .ne. 0) then
            print *, 'ERROR: Failed to open ', file, ': iostat = ', iostat
            print *, 'Exiting...'
            stop
        end if
    end subroutine check_file_open
    subroutine open_file(unit, form, file, status, action)
        implicit none
        character(len=*), intent(in) :: form, file, status, action
        integer, intent(inout) :: unit
        character(:), allocatable :: file_status
        integer :: iostat
        call search_unused_file_unit(unit)
        file_status = trim(status)
        call lowercase(file_status)
        if (file_status /= 'old' .and. file_status /= 'new' .and. file_status /= 'replace') then
            print *, 'ERROR: file_status must be old, new or replace. file_status = ', file_status
            print *, 'Exiting...'
            stop
        end if
        open (unit, form=form, file=file, status=status, iostat=iostat, action=action)
        call check_file_open(file, iostat)
    end subroutine open_file

    subroutine open_unformatted_file(unit, file, status, optional_action)
        implicit none
        character(len=*), intent(in), optional :: optional_action
        character(len=*), intent(in) :: file, status
        integer, intent(inout) :: unit
        character(:), allocatable :: actual_action, trimmed_action, form

        if (present(optional_action)) then
            trimmed_action = trim(optional_action)
            call lowercase(trimmed_action)
            if (trimmed_action /= 'read' .and. trimmed_action /= 'write' .and. trimmed_action /= 'readwrite') then
                print *, 'ERROR: trimmed_action must be read, write or readwrite. trimmed_action = ', trimmed_action
                print *, 'FILE NAME: ', file
                print *, 'Exiting...'
                stop
            end if
            actual_action = trimmed_action
        else
            actual_action = 'readwrite'
        end if
        form = 'unformatted'
        call open_file(unit=unit, form=form, file=file, status=status, action=actual_action)
    end subroutine open_unformatted_file

    subroutine open_fomatted_file(unit, file, status, optional_action)
        implicit none
        character(len=*), intent(in), optional :: optional_action
        character(len=*), intent(in) :: file, status
        integer, intent(inout) :: unit
        character(:), allocatable :: form, actual_action, trimmed_action

        if (present(optional_action)) then
            trimmed_action = trim(optional_action)
            call lowercase(trimmed_action)
            if (trimmed_action /= 'read' .and. trimmed_action /= 'write' .and. trimmed_action /= 'readwrite') then
                print *, 'ERROR: trimmed_action must be read, write or readwrite. trimmed_action = ', trimmed_action
                print *, 'FILE NAME: ', file
                print *, 'Exiting...'
                stop
            end if
            actual_action = trimmed_action
        else
            actual_action = 'readwrite'
        end if
        form = 'formatted'
        call open_file(unit=unit, form=form, file=file, status=status, action=actual_action)
    end subroutine open_fomatted_file
end module module_file_manager
