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
    print *, "open"
    open (2, file='result', form="formatted")
    print *, 'before write'
    write (2, *) string
    print *, 'end write'
    close (2)
end program main
