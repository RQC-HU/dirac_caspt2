module file_manage_module
    implicit none
!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
! file_manage_module
! Copyright (c) by the authors of rel-caspt2.
! Author K.Noda
!
! This is a utility module that manages the file unit number.
!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
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

    subroutine open_file(file_name, file_id)
        implicit none
        character(len=*), intent(in) :: file_name
        integer, intent(inout) :: file_id

        call search_unused_file_unit(file_id)
        open (file_id, file=file_name, status='unknown')
    end subroutine open_file
end module file_manage_module
