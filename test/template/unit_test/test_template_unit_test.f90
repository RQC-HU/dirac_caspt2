program main
    use module_file_manager
    implicit none
    character(:), allocatable :: string
    integer :: unit_new

    string = "Hello World"
    call open_formatted_file(unit=unit_new, file='result.out', status="replace", optional_action='write')
    write (unit_new, *) string
    close (unit_new)
end program main
