! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE tramo1(i, j, int1)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables
    use module_index_utils, only: get_mo_range
    use module_realonly, only: realonly

    Implicit NONE
    integer, intent(in)  :: i, j
    complex*16, intent(out) :: int1

    integer :: i0, j0, sym1, sym2, count_n_list1, count_n_list2, idx_i0, idx_j0
    integer :: start_idx_i, start_idx_j, end_idx_i, end_idx_j
    integer, allocatable :: n_list_sym1(:), n_list_sym2(:)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    int1 = 0.0d+00
    sym1 = irpamo(i)
    sym2 = irpamo(j)

    If (MULTB_D(sym1, sym2) == 1) then

        call get_mo_range(i, start_idx_i, end_idx_i)
        call get_mo_range(j, start_idx_j, end_idx_j)

        count_n_list1 = count(irpamo(start_idx_i:end_idx_i) == sym1)
        allocate (n_list_sym1(count_n_list1))
        call craete_n_list(start_idx_i, end_idx_i, sym1, n_list_sym1)

        count_n_list2 = count(irpamo(start_idx_j:end_idx_j) == sym2)
        allocate (n_list_sym2(count_n_list2))
        call craete_n_list(start_idx_j, end_idx_j, sym2, n_list_sym2)

        do idx_j0 = 1, count_n_list2
            do idx_i0 = 1, count_n_list1
                i0 = n_list_sym1(idx_i0)
                j0 = n_list_sym2(idx_j0)
                if (realonly%is_realonly()) then
                    int1 = int1 + fock_real(i0, i)*one_elec_int_r(i0, j0)*fock_real(j0, j)
                else
                    int1 = int1 + DCONJG(fock_cmplx(i0, i))*DCMPLX(one_elec_int_r(i0, j0), one_elec_int_i(i0, j0))*fock_cmplx(j0, j)
                end if
            end do
        end do

    End if
contains
    subroutine craete_n_list(idx_start, idx_end, isym, n_list)
        implicit none
        integer, intent(in) :: idx_start, idx_end, isym
        integer, intent(out) :: n_list(:)
        integer :: idx, idx_list(idx_start:idx_end)

        do idx = idx_start, idx_end
            idx_list(idx) = idx
        end do
        n_list = pack(idx_list, mask=irpamo(idx_start:idx_end) == isym)

    end subroutine craete_n_list

End subroutine tramo1

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE tramo2(i, j, k, l, int2)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables
    use module_index_utils, only: get_mo_range
    use module_realonly, only: realonly

    Implicit NONE
    integer, intent(in) :: i, j, k, l

    complex*16, intent(out) :: int2

    integer     :: i0, j0, k0, l0, sym1, sym2, sym3, sym4, sym5, sym6
    integer     :: idx_i0, idx_j0, idx_k0, idx_l0, count_sym1_i, count_sym2_j, count_sym3_k, count_sym4_l
    integer     :: start_idx_i, start_idx_j, start_idx_k, start_idx_l, end_idx_i, end_idx_j, end_idx_k, end_idx_l
    integer, allocatable    :: list_sym1_i(:), list_sym2_j(:), list_sym3_k(:), list_sym4_l(:)
    real(8)      :: i2r, i2i
    complex*16  :: cmplxint

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    int2 = 0.0d+00
    cmplxint = 0.0d+00
    sym1 = irpamo(i)
    sym2 = irpamo(j)
    sym3 = irpamo(k)
    sym4 = irpamo(l)

    sym5 = MULTB_D(sym1, sym2)
    sym6 = MULTB_D(sym3, sym4)

    If (MULTB_S(sym5, sym6) == 1) then

        call get_mo_range(i, start_idx_i, end_idx_i)
        call get_mo_range(j, start_idx_j, end_idx_j)
        call get_mo_range(k, start_idx_k, end_idx_k)
        call get_mo_range(l, start_idx_l, end_idx_l)

        count_sym1_i = count(irpamo(start_idx_i:end_idx_i) == sym1)
        allocate (list_sym1_i(count_sym1_i))
        call craete_n_list(start_idx_i, end_idx_i, sym1, list_sym1_i)

        count_sym2_j = count(irpamo(start_idx_j:end_idx_j) == sym2)
        allocate (list_sym2_j(count_sym2_j))
        call craete_n_list(start_idx_j, end_idx_j, sym2, list_sym2_j)

        count_sym3_k = count(irpamo(start_idx_k:end_idx_k) == sym3)
        allocate (list_sym3_k(count_sym3_k))
        call craete_n_list(start_idx_k, end_idx_k, sym3, list_sym3_k)

        count_sym4_l = count(irpamo(start_idx_l:end_idx_l) == sym4)
        allocate (list_sym4_l(count_sym4_l))
        call craete_n_list(start_idx_l, end_idx_l, sym4, list_sym4_l)

        do idx_l0 = 1, count_sym4_l
            do idx_k0 = 1, count_sym3_k
                do idx_j0 = 1, count_sym2_j
                    do idx_i0 = 1, count_sym1_i
                        i0 = list_sym1_i(idx_i0); j0 = list_sym2_j(idx_j0); k0 = list_sym3_k(idx_k0); l0 = list_sym4_l(idx_l0)

                        if (realonly%is_realonly()) then
                            i2r = inttwr(i0, j0, k0, l0)
                            cmplxint = DCMPLX(i2r, 0.0d+00)
                            int2 = int2 + fock_real(i0, i)*fock_real(k0, k)*fock_real(j0, j)*fock_real(l0, l)*cmplxint
                        else
                            i2r = inttwr(i0, j0, k0, l0)
                            i2i = inttwi(i0, j0, k0, l0)
                            cmplxint = DCMPLX(i2r, i2i)

                            int2 = int2 + DCONJG(fock_cmplx(i0, i))*DCONJG(fock_cmplx(k0, k))* &
                                   fock_cmplx(j0, j)*fock_cmplx(l0, l)*cmplxint
                        end if
                    end do
                end do
            end do
        end do

        deallocate (list_sym1_i, list_sym2_j, list_sym3_k, list_sym4_l)

    End if
contains
    subroutine craete_n_list(idx_start, idx_end, isym, n_list)
        implicit none
        integer, intent(in) :: idx_start, idx_end, isym
        integer, intent(out) :: n_list(:)
        integer :: idx, idx_list(idx_start:idx_end)

        do idx = idx_start, idx_end
            idx_list(idx) = idx
        end do
        ! n_list consists of the list that satisfies the mask condition
        n_list = pack(idx_list, mask=irpamo(idx_start:idx_end) == isym)

    end subroutine craete_n_list

End subroutine tramo2
