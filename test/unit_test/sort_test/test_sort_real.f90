program main
    use module_file_manager
    use module_sort_swap
    implicit none
    real(8) :: want_to_sort_real(6) = (/8.1, -9.2, 10000.58, -897.0, 123456789.0, 0.0000000010/)
    integer :: unit_new
    call heapSort(want_to_sort_real, .false.)
    call open_formatted_file(unit=unit_new, file="real.out",status='replace' ,optional_action="write")
    print *, want_to_sort_real
    write (unit_new, *) want_to_sort_real
    close (unit_new)
end program main
