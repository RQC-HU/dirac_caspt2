program main
    use four_caspt2_module, only: ras3_list
    use read_input_module
    implicit none
    open (5, file='input', form='formatted')
    call ras_read(ras3_list, 3)
    close (5)
    open (2, file='result.out', form="formatted")
    write (2, *) ras3_list
    close (2)
end program main
