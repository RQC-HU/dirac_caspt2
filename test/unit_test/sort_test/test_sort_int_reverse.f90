program main
    use module_file_manager
    use module_sort_swap
    implicit none
    integer :: want_to_sort(21) = (/8, 9, 10, 11, 12, 13, 14, 15, 16, 169, 170, 171, 172, 173, 174, 175, 1, 3, 5, 156, 189/)
    integer :: new_unit = 20
    call heapSort(want_to_sort, .true.)
    call open_formatted_file(unit=new_unit, file="int_reverse.out",status='replace' ,optional_action="write")
    write (new_unit, *) want_to_sort
    close (new_unit)
end program main
