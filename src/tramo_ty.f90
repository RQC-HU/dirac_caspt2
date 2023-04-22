! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE tramo1_ty(i, j, int1)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

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
                int1 = int1 + DCONJG(f(i0, i))*DCMPLX(one_elec_int_r(i0, j0), one_elec_int_i(i0, j0))*f(j0, j)
            end do
        end do

    End if
contains
    subroutine get_mo_range(mo_idx, start_idx, end_idx)
        implicit none
        integer, intent(in) :: mo_idx
        integer, intent(out) :: start_idx, end_idx
        if (mo_idx <= ninact) then
            start_idx = 1
            end_idx = ninact
        elseif (mo_idx <= ninact + nact) then
            start_idx = ninact + 1
            end_idx = ninact + nact
        elseif (mo_idx <= ninact + nact + nsec) then
            start_idx = ninact + nact + 1
            end_idx = ninact + nact + nsec
        else
            print *, "invalid mo_idx = ", mo_idx
        end if

    end subroutine get_mo_range
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

End subroutine tramo1_ty

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE tramo2_ty(i, j, k, l, int2)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module
    use module_error, only: stop_with_errorcode

    Implicit NONE
    integer, intent(in) :: i, j, k, l

    complex*16, intent(out) :: int2

    integer     :: i0, j0, k0, l0, sym1, sym2, sym3, sym4, sym5, sym6
    integer     :: idx_i0, idx_j0, idx_k0, idx_l0, count_sym1_i, count_sym2_j, count_sym3_k, count_sym4_l
    integer     :: start_idx_i, start_idx_j, start_idx_k, start_idx_l, end_idx_i, end_idx_j, end_idx_k, end_idx_l
    integer, allocatable    :: list_sym1_i(:), list_sym2_j(:), list_sym3_k(:), list_sym4_l(:)
    real*8      :: i2r, i2i
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

                        i2r = inttwr(i0, j0, k0, l0)
                        i2i = inttwi(i0, j0, k0, l0)
                        cmplxint = DCMPLX(i2r, i2i)

                        int2 = int2 + DCONJG(f(i0, i))*DCONJG(f(k0, k))*f(j0, j)*f(l0, l)*cmplxint

                    end do
                end do
            end do
        end do

        deallocate (list_sym1_i, list_sym2_j, list_sym3_k, list_sym4_l)

    End if
contains
    subroutine get_mo_range(mo_idx, mo_start, mo_end)
        implicit none
        integer, intent(in) :: mo_idx
        integer, intent(out) :: mo_start, mo_end

        if (mo_idx <= ninact) then
            mo_start = 1
            mo_end = ninact
        elseif (mo_idx <= ninact + nact) then
            mo_start = ninact + 1
            mo_end = ninact + nact
        elseif (mo_idx <= ninact + nact + nsec) then
            mo_start = ninact + nact + 1
            mo_end = ninact + nact + nsec
        else
            print *, "invalid mo_idx =", mo_idx
            call stop_with_errorcode(1)
        end if

    end subroutine get_mo_range
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

End subroutine tramo2_ty
