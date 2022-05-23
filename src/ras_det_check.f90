module ras_det_check
    implicit none
    private
    public ras1_det_check
contains
    function ras1_det_check(i,upper_allowed_hole) result(is_det_allowed)
        ! function ras1_det_check(i,upper_allowed_hole) result(is_det_allowed)
        ! This function returns true if the determinant (i) is allowed
        integer, intent(in) :: i,upper_allowed_hole
        integer :: num_of_electron
        logical :: is_det_allowed
        call conunt_num_of_elec(i, num_of_electron)
        print *, is_det_allowed
        is_det_allowed = 2 - upper_allowed_hole <= num_of_electron  .and. num_of_electron <= 2
        print *, is_det_allowed
    end function ras1_det_check
    function ras3_det_check(i,upper_allowed_electron) result(is_det_allowed)
        ! function ras3_det_check(i,upper_allowed_electron) result(is_det_allowed)
        ! This function returns true if the determinant (i) is allowed
        integer, intent(in) :: i,upper_allowed_electron
        integer :: num_of_electron
        logical :: is_det_allowed
        call conunt_num_of_elec(i, num_of_electron)
        print *, is_det_allowed
        is_det_allowed = num_of_electron <= upper_allowed_electron
        print *, is_det_allowed
    end function ras3_det_check

    subroutine conunt_num_of_elec(i, num_of_electron)
        use four_caspt2_module, only: spinor_num_ras1

        implicit none
        integer, intent(in) ::  i
        integer, intent(out) :: num_of_electron
        print *, 2**spinor_num_ras1-1, i, num_of_electron
        num_of_electron = ras_bit_calculate(i)
        print *, 2**spinor_num_ras1-1, i, num_of_electron


    end subroutine conunt_num_of_elec

    function ras_bit_calculate(determinant) result(num_of_electron)
        use four_caspt2_module, only: spinor_num_ras1
        implicit none
        integer, intent(in) :: determinant
        integer :: num_of_electron, multiply, ras1_bit
        ras1_bit = 2 ** spinor_num_ras1 - 1
        ! ras_bitとdeterminantとのbit論理積
        multiply = iand(ras1_bit, determinant)
        num_of_electron = popcnt(multiply)

    end function ras_bit_calculate
end module ras_det_check
