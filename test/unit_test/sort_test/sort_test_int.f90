program main
    use module_sort_swap
    implicit none
    integer :: want_to_sort(21) = (/8, 9, 10, 11, 12, 13, 14, 15, 16, 169, 170, 171, 172, 173, 174, 175, 1, 3, 5, 156, 189/)
    call heapSort(want_to_sort, .false.)
    open (1, file='intout', form='formatted')
    write (1, *) want_to_sort
    close (1)
end program main
