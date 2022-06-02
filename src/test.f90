program test
    use ras_check_det, only: ras1_det_check

    implicit none
    integer :: i,ras_bit,upper_allowed_hole,ras1_bit
    logical :: is_allow
    i = 0
    ras_bit = 3
    upper_allowed_hole = 1
    print *, ras1_bit,i,upper_allowed_hole
    is_allow = ras1_det_check(i,ras_bit,upper_allowed_hole)

end program test
