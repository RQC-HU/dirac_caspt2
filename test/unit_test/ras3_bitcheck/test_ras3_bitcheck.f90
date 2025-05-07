program ras3_bitcheck
    use, intrinsic :: iso_fortran_env, only: int64
    use module_global_variables
    use module_file_manager
    use read_input_module
    use ras_det_check
    implicit none
    integer :: i, unit_new
    logical :: is_allow

    call open_formatted_file(unit=unit_new, file='active.inp', status="old", optional_action='read')
    call read_input(unit_new, .true.)
    close (unit_new)

    call open_formatted_file(unit=unit_new, file='result.out', status="replace", optional_action='write')
    do i = 1, 2**nact - 1
        is_allow = satisfy_ras3_condition(int(i, kind=int64), int(ras3_max_elec, kind=int64))
        if (is_allow) then
            print '(i4,b20)', i, i
            write (unit_new, '(i4,b20)') i, i
        end if
    end do
    close (unit_new)
end program ras3_bitcheck
