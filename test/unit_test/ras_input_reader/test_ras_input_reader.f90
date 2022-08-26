program main
    use four_caspt2_module, only: ras3_list
    use module_file_manager
    use read_input_module
    implicit none
    integer :: new_unit = 20
    call open_formatted_file(unit=new_unit, file='input', status='old', optional_action='read')
    call ras_read(new_unit, ras3_list, 3)
    close (new_unit)
    call open_formatted_file(unit=new_unit, file='result.out', status='replace', optional_action='write')
    write (new_unit, *) ras3_list
    close (new_unit)
end program main
