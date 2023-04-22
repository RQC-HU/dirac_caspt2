program main
    use module_file_manager
    use module_sort_swap
    implicit none
    integer :: want_to_sort(21) = (/8, 9, 10, 11, 12, 13, 14, 15, 16, 169, 170, 171, 172, 173, 174, 175, 1, 3, 5, 156, 189/)
    integer :: unit_new
    call heapSort(want_to_sort, .false.)
    call open_formatted_file(unit=unit_new, file="int.out",status='replace' ,optional_action="write")
    write (unit_new, *) want_to_sort
    close (unit_new)
end program main
