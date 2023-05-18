program main
    use module_file_manager
    use read_input_module
    implicit none
    character(len=100) :: input
    character(:), allocatable :: string
    integer :: unit_new
    call open_formatted_file(unit=unit_new, file="input", status='old', optional_action='read')
    read (unit_new, '(a)') input
    close (unit_new)

    allocate (string, source=trim(input))
    call uppercase(string)
    call open_formatted_file(unit=unit_new, file="result.out", status='replace', optional_action='write')
    write (unit_new, *) string
    close (unit_new)
end program main
