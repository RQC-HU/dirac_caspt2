program main
    use read_input_module
    implicit none
    character(100) :: input
    integer :: count
    count = 1
    open (5, file='input', form='formatted')
    read (5,'(a)') input
    print *, input
    close (5)
    call lowercase(input)
    print *, "open"
    open (2, file='output', form="formatted")
    print *, 'before write'
    write (2, *) input
    print *, 'end write'
    close (2)
end program main
