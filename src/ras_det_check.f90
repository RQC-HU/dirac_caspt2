module ras_det_check
    use four_caspt2_module, only: rank, ras1_list, ras3_list
    implicit none
    private
    public ras1_det_check, ras3_det_check
contains
    function ras1_det_check(i, upper_allowed_hole) result(is_det_allowed)
        ! function ras1_det_check(i,upper_allowed_hole) result(is_det_allowed)
        ! This function returns true if the determinant (i) is allowed
        integer, intent(in) :: i, upper_allowed_hole
        integer :: num_of_electron, ras1_bit
        logical :: is_det_allowed
        ras1_bit = 2**size(ras1_list, 1) - 1
        call conunt_num_of_elec(i, ras1_bit, num_of_electron)
        is_det_allowed = 2 - upper_allowed_hole <= num_of_electron .and. num_of_electron <= 1
    end function ras1_det_check
    function ras3_det_check(i, upper_allowed_electron) result(is_det_allowed)
        ! function ras3_det_check(i,upper_allowed_electron) result(is_det_allowed)
        ! This function returns true if the determinant (i) is allowed
        integer, intent(in) :: i, upper_allowed_electron
        integer :: num_of_electron, ras3_bit
        logical :: is_det_allowed
        ras3_bit = ishft(2**size(ras3_list, 1) - 1, size(ras1_list, 1))
        call conunt_num_of_elec(i, ras3_bit, num_of_electron)
        ! print *, 'res', i, num_of_electron
        is_det_allowed = num_of_electron <= upper_allowed_electron
    end function ras3_det_check

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
