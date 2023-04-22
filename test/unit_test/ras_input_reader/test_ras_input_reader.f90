program main
    use four_caspt2_module, only: ras3_list
    use module_file_manager
    use read_input_module
    implicit none
    integer :: unit_new
    call open_formatted_file(unit=unit_new, file='input', status='old', optional_action='read')
    call ras_read(unit_new, ras3_list, 3)
    close (unit_new)
    call open_formatted_file(unit=unit_new, file='result.out', status='replace', optional_action='write')
    write (unit_new, *) ras3_list
    close (unit_new)
end program main
