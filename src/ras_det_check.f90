module ras_det_check
    implicit none
    private
    public ras1_det_check
contains
    function ras1_det_check(i,ras_bit,upper_allowed_hole) result(is_det_allowed)
        ! function ras1_det_check(i,ras_bit,upper_allowed_hole) result(is_det_allowed)
        ! This function returns true if the determinant (i) is allowed
        integer, intent(in) :: i,ras_bit,upper_allowed_hole
        integer :: num_of_electron
        logical :: is_det_allowed
        call conunt_num_of_elec(i, ras_bit, upper_allowed_hole, num_of_electron)
        print *, is_det_allowed
        is_det_allowed = 2 - upper_allowed_hole <= num_of_electron  .and. num_of_electron <= 2
        print *, is_det_allowed
    end function ras1_det_check
    function ras3_det_check(i,ras_bit,upper_allowed_electron) result(is_det_allowed)
        ! function ras3_det_check(i,ras_bit,upper_allowed_electron) result(is_det_allowed)
        ! This function returns true if the determinant (i) is allowed
        integer, intent(in) :: i,ras_bit,upper_allowed_electron
        integer :: num_of_electron
        logical :: is_det_allowed
        call conunt_num_of_elec(i, ras_bit, upper_allowed_electron, num_of_electron)
        print *, is_det_allowed
        is_det_allowed = num_of_electron <= upper_allowed_electron
        print *, is_det_allowed
    end function ras3_det_check

    subroutine conunt_num_of_elec(i, ras_bit, upper_allowed_hole, num_of_electron)

        implicit none
        integer, intent(in) ::  i, ras_bit, upper_allowed_hole
        integer, intent(out) :: num_of_electron
        print *, ras_bit, i, upper_allowed_hole, num_of_electron
        num_of_electron = ras_bit_calculate(ras_bit, i)
        print *, ras_bit, i, upper_allowed_hole, num_of_electron


    end subroutine conunt_num_of_elec

    function ras_bit_calculate(ras_bit, determinant) result(num_of_electron)
        implicit none
        integer, intent(in) :: ras_bit, determinant
        integer :: num_of_electron, multiply

        ! ras_bitとdeterminantとのbit論理積
        multiply = iand(ras_bit, determinant)
        num_of_electron = popcnt(multiply)

    end function ras_bit_calculate
end module ras_det_check
