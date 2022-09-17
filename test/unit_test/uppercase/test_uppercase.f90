program main
    use module_file_manager
    use read_input_module
    implicit none
    character(100) :: input
    character(:), allocatable :: string
    integer :: new_unit = 20
    call open_formatted_file(unit=new_unit, file="input", status='old', optional_action='read')
    read (new_unit, '(a)') input
    close (new_unit)

    allocate (string, source=trim(input))
    call uppercase(string)
    call open_formatted_file(unit=new_unit, file="result.out", status='replace', optional_action='write')
    write (new_unit, *) string
    close (new_unit)
end program main
