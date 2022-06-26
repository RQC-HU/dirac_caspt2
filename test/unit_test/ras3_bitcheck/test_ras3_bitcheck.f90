program ras3_bitcheck
    use four_caspt2_module
    use read_input_module
    use ras_det_check
    implicit none
    integer :: i, shift
    logical :: is_allow

    call read_input
    open (10, file="result", form="formatted")
    do i = 1, 2**nact - 1
        shift = ishft(i, -ras1num)
        is_allow = ras3_det_check(shift, ras3_max_elec)
        if (is_allow) then
            print '(i4,b20)', i, i
            write (10, '(i4,b20)'), i, i
        end if
    end do
    close (10)
end program ras3_bitcheck
