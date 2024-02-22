program main
    use module_file_manager
    use module_time
    implicit none
    integer :: unit_new, num_times, i

    call open_formatted_file(unit_new, "input", "old", "read")

    read (unit_new, *) num_times

    do i = 1, num_times
        read (unit_new, *) start_time%year, start_time%month, start_time%date, &
            start_time%hour, start_time%min, start_time%sec, start_time%msec
        read (unit_new, *) end_time%year, end_time%month, end_time%date, &
            end_time%hour, end_time%min, end_time%sec, end_time%msec
        call print_time_diff(start_time, end_time)
    end do

end program main
