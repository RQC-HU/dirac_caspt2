program main
    use module_sort_swap
    implicit none
    real(8) :: want_to_sort_real(6) = (/8.1, -9.2, 10000.58, -897.0, 123456789.0, 0.0000000010/)
    call heapSort(want_to_sort_real, .false.)
    open (1, file='realout', form='formatted')
    print *, want_to_sort_real
    write (1, *) want_to_sort_real
    close (1)
end program main
