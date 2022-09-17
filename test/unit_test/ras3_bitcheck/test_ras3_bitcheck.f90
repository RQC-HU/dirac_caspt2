program ras3_bitcheck
    use four_caspt2_module
    use module_file_manager
    use read_input_module
    use ras_det_check
    implicit none
    integer :: i, new_unit
    logical :: is_allow

    new_unit = 20
    call open_formatted_file(unit=new_unit, file='active.inp', status="old", optional_action='read')
    call read_input(new_unit)
    close (new_unit)

    call open_formatted_file(unit=new_unit, file='result', status="replace", optional_action='write')
    do i = 1, 2**nact - 1
        is_allow = ras3_det_check(i, ras3_max_elec)
        if (is_allow) then
            print '(i4,b20)', i, i
            write (new_unit, '(i4,b20)') i, i
        end if
    end do
    close (new_unit)
end program ras3_bitcheck
