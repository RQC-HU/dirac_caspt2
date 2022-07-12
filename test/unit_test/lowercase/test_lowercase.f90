program main
    use read_input_module
    implicit none
    character(100) :: input
    character(:), allocatable :: string
    integer :: count
    count = 1
    open (5, file='input', form='formatted')
    read (5, '(a)') input
    string = trim(input)
    close (5)
    call lowercase(string)
    open (2, file='result.out', form="formatted")
    write (2, *) string
    close (2)
end program main
