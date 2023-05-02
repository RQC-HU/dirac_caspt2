module ras_det_check
    use four_caspt2_module, only: rank, ras1_list, ras2_list
    implicit none
    private
    public satisfy_ras1_condition, satisfy_ras3_condition
contains
    function satisfy_ras1_condition(i, upper_allowed_hole) result(is_det_allowed)
        ! function satisfy_ras1_condition(i,upper_allowed_hole) result(is_det_allowed)
        ! This function returns true if the determinant (i) is allowed
        use four_caspt2_module, only: ras1_size, min_hole_ras1
        integer, intent(in) :: i, upper_allowed_hole
        integer :: num_of_electron, ras1_bit
        logical :: is_det_allowed
        ras1_bit = 2**ras1_size - 1
        call conunt_num_of_elec(i, ras1_bit, num_of_electron)
        is_det_allowed = ras1_size - upper_allowed_hole <= num_of_electron .and. num_of_electron <= ras1_size - min_hole_ras1
    end function satisfy_ras1_condition
    function satisfy_ras3_condition(i, upper_allowed_electron) result(is_det_allowed)
        ! function satisfy_ras3_condition(i,upper_allowed_electron) result(is_det_allowed)
        ! This function returns true if the determinant (i) is allowed
        use four_caspt2_module, only: ras1_size, ras2_size
        integer, intent(in) :: i, upper_allowed_electron
        integer :: num_of_electron, ras3_bit, width_of_shift
        logical :: is_det_allowed
        ras3_bit = i
        width_of_shift = 0
        if (ras1_size /= 0) then
            ras3_bit = ishft(ras3_bit, -ras1_size)
            width_of_shift = width_of_shift + ras1_size
        end if
        if (ras2_size /= 0) then
            ras3_bit = ishft(ras3_bit, -ras2_size)
            width_of_shift = width_of_shift + ras2_size
        end if
        ras3_bit = ishft(ras3_bit, width_of_shift)
        call conunt_num_of_elec(i, ras3_bit, num_of_electron)
        is_det_allowed = num_of_electron <= upper_allowed_electron
    end function satisfy_ras3_condition

    subroutine conunt_num_of_elec(i, bit, num_of_electron)
        implicit none
        integer, intent(in) ::  i, bit
        integer, intent(out) :: num_of_electron
        num_of_electron = ras_bit_calculate(i, bit)

    end subroutine conunt_num_of_elec

    function ras_bit_calculate(determinant, bit) result(num_of_electron)
        implicit none
        integer, intent(in) :: determinant, bit
        integer :: num_of_electron, multiply

        ! ras_bitとdeterminantとのbit論理積
        multiply = iand(bit, determinant)
        num_of_electron = popcnt(multiply)

    end function ras_bit_calculate
end module ras_det_check
