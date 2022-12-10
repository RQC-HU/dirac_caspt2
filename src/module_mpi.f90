module module_mpi
    use four_caspt2_module, only: rank, nprocs, ierr, max_i4
    implicit none
#ifdef HAVE_MPI
    include 'mpif.h'
    private
    public allreduce_wrapper
    interface allreduce_wrapper
        ! Default operation is MPI_SUM
        module procedure allreduce_i, allreduce_i_1, allreduce_i_2, allreduce_i_3, allreduce_i_4, &
            allreduce_r, allreduce_r_1, allreduce_r_2, allreduce_r_3, allreduce_r_4, &
            allreduce_c, allreduce_c_1, allreduce_c_2, allreduce_c_3, allreduce_c_4
    end interface allreduce_wrapper
    integer, parameter, public :: op_mpi_max = MPI_MAX          ! MPI_MAX      最大値
    integer, parameter, public :: op_mpi_min = MPI_MIN          ! MPI_MIN      最小値
    integer, parameter, public :: op_mpi_sum = MPI_SUM          ! MPI_SUM      和
    integer, parameter, public :: op_mpi_prod = MPI_PROD        ! MPI_PROD     積
    integer, parameter, public :: op_mpi_land = MPI_LAND        ! MPI_LAND     論理積
    integer, parameter, public :: op_mpi_band = MPI_BAND        ! MPI_BAND     ビット演算の積
    integer, parameter, public :: op_mpi_lor = MPI_LOR          ! MPI_LOR      論理和
    integer, parameter, public :: op_mpi_bor = MPI_BOR          ! MPI_BOR      ビット演算の和
    integer, parameter, public :: op_mpi_lxor = MPI_LXOR        ! MPI_LXOR     排他的論理和
    integer, parameter, public :: op_mpi_bxor = MPI_BXOR        ! MPI_BXOR     ビット演算の排他的論理和
    integer, parameter, public :: op_mpi_maxloc = MPI_MAXLOC    ! MPI_MAXLOC   最大値と位置
    integer, parameter, public :: op_mpi_minloc = MPI_MINLOC    ! MPI_MINLOC   最小値と位置

    integer, parameter :: ops(12) = (/op_mpi_max, op_mpi_min, op_mpi_sum, op_mpi_prod, &
                                      op_mpi_land, op_mpi_band, op_mpi_lor, op_mpi_bor, &
                                      op_mpi_lxor, op_mpi_bxor, op_mpi_maxloc, op_mpi_minloc/)
contains

    subroutine allreduce_i(mat, optional_op)
        ! Allreduce for a integer value
        implicit none
        integer, intent(inout) :: mat
        integer, optional, intent(in) :: optional_op
        integer :: ii, ie
        integer :: i, idx_end, cnt
        integer :: op = op_mpi_sum ! default operation
        integer :: datatype

        if (present(optional_op)) then
            call check_operation(optional_op)
            op = optional_op
        end if
        if (sizeof(mat) == 4) then
            datatype = MPI_INTEGER4 ! 4 byte integer
        else
            datatype = MPI_INTEGER8 ! 8 byte integer
        end if

        call MPI_Allreduce(MPI_IN_PLACE, mat, 1, datatype, op, MPI_COMM_WORLD, ierr)
        call check_ierr(ierr)

    end subroutine allreduce_i

    subroutine allreduce_i_1(mat, optional_op)
        ! Allreduce for 1 dimensional integer array
        implicit none
        integer, intent(inout) :: mat(:)
        integer, optional, intent(in) :: optional_op
        integer :: ii, ie
        integer :: i, idx_end, cnt
        integer :: op = op_mpi_sum ! default operation
        integer :: datatype

        if (present(optional_op)) then
            call check_operation(optional_op)
            op = optional_op
        end if

        ! Set the first and last index of the array for each dimension.
        ii = lbound(mat, 1); ie = ubound(mat, 1)
        if (sizeof(mat(ii)) == 4) then
            datatype = MPI_INTEGER4 ! 4 byte integer
        else
            datatype = MPI_INTEGER8 ! 8 byte integer
        end if

        ! Because of the limitation of the size of the array, the array is divided into several parts and allreduce is performed.
        if (max_i4 < size(mat, 1)) then ! 1st index is larger than the limit(max_i4)
            do i = ii, ie, max_i4
                idx_end = min(i + max_i4 - 1, ie)
                cnt = idx_end - i + 1
                call MPI_Allreduce(MPI_IN_PLACE, mat(i:idx_end), &
                                   cnt, datatype, op, MPI_COMM_WORLD, ierr)
                call check_ierr(ierr)
            end do
        else ! no limit
            call MPI_Allreduce(MPI_IN_PLACE, mat, size(mat), datatype, op, MPI_COMM_WORLD, ierr)
            call check_ierr(ierr)
        end if
    end subroutine allreduce_i_1

    subroutine allreduce_i_2(mat, optional_op)
        ! Allreduce for 2 dimensional integer array
        implicit none
        integer, intent(inout) :: mat(:, :)
        integer, optional, intent(in) :: optional_op
        integer :: ii, ie, ji, je
        integer :: i, j, idx_end, cnt
        integer :: op = op_mpi_sum ! default operation
        integer :: datatype

        if (present(optional_op)) then
            call check_operation(optional_op)
            op = optional_op
        end if

        ! Set the first and last index of the array for each dimension.
        ii = lbound(mat, 1); ie = ubound(mat, 1)
        ji = lbound(mat, 2); je = ubound(mat, 2)
        if (sizeof(mat(ii, ji)) == 4) then
            datatype = MPI_INTEGER4 ! 4 byte integer
        else
            datatype = MPI_INTEGER8 ! 8 byte integer
        end if
        ! Because of the limitation of the size of the array, the array is divided into several parts and allreduce is performed.
        if (max_i4 < size(mat, 1)) then ! 1st index is larger than the limit(max_i4)
            do j = ji, je
                do i = ii, ie, max_i4
                    idx_end = min(i + max_i4 - 1, ie)
                    cnt = idx_end - i + 1
                    call MPI_Allreduce(MPI_IN_PLACE, mat(i:idx_end, j), &
                                       cnt, datatype, op, MPI_COMM_WORLD, ierr)
                    call check_ierr(ierr)
                end do
            end do
        else if (max_i4 < size(mat, 1)*size(mat, 2)) then ! 1st index * 2nd index is larger than the limit(max_i4)
            do j = ji, je, max_i4
                idx_end = min(j + max_i4 - 1, je)
                cnt = (idx_end - j + 1)*size(mat, 1)
                call MPI_Allreduce(MPI_IN_PLACE, mat(:, j:idx_end), &
                                   cnt, datatype, op, MPI_COMM_WORLD, ierr)
                call check_ierr(ierr)
            end do
        else ! no limit
            call MPI_Allreduce(MPI_IN_PLACE, mat, size(mat), datatype, op, MPI_COMM_WORLD, ierr)
            call check_ierr(ierr)
        end if
    end subroutine allreduce_i_2

    subroutine allreduce_i_3(mat, optional_op)
        ! Allreduce for 3 dimensional integer array
        implicit none
        integer, intent(inout) :: mat(:, :, :)
        integer, optional, intent(in) :: optional_op
        integer :: ii, ie, ji, je, ki, ke
        integer :: i, j, k, idx_end, cnt
        integer :: op = op_mpi_sum ! default operation
        integer :: datatype

        if (present(optional_op)) then
            call check_operation(optional_op)
            op = optional_op
        end if

        ! Set the first and last index of the array for each dimension.
        ii = lbound(mat, 1); ie = ubound(mat, 1)
        ji = lbound(mat, 2); je = ubound(mat, 2)
        ki = lbound(mat, 3); ke = ubound(mat, 3)
        if (sizeof(mat(ii, ji, ki)) == 4) then
            datatype = MPI_INTEGER4 ! 4 byte integer
        else
            datatype = MPI_INTEGER8 ! 8 byte integer
        end if
        ! Because of the limitation of the size of the array, the array is divided into several parts and allreduce is performed.
        if (max_i4 < size(mat, 1)) then ! 1st index is larger than the limit(max_i4)
            do k = ki, ke
                do j = ji, je
                    do i = ii, ie, max_i4
                        idx_end = min(i + max_i4 - 1, ie)
                        cnt = idx_end - i + 1
                        call MPI_Allreduce(MPI_IN_PLACE, mat(i:idx_end, j, k), &
                                           cnt, datatype, op, MPI_COMM_WORLD, ierr)
                        call check_ierr(ierr)
                    end do
                end do
            end do
        else if (max_i4 < size(mat, 1)*size(mat, 2)) then ! 1st index * 2nd index is larger than the limit(max_i4)
            do k = ki, ke
                do j = ji, je, max_i4
                    idx_end = min(j + max_i4 - 1, je)
                    cnt = (idx_end - j + 1)*size(mat, 1)
                    call MPI_Allreduce(MPI_IN_PLACE, mat(:, j:idx_end, k), &
                                       cnt, datatype, op, MPI_COMM_WORLD, ierr)
                    call check_ierr(ierr)
                end do
            end do
        else if (max_i4 < size(mat, 1)*size(mat, 2)*size(mat, 3)) then ! 1st index * 2nd index * 3rd index is larger than the limit(max_i4)
            do k = ki, ke, max_i4
                idx_end = min(k + max_i4 - 1, ke)
                cnt = (idx_end - k + 1)*size(mat, 1)*size(mat, 2)
                call MPI_Allreduce(MPI_IN_PLACE, mat(:, :, k:idx_end), &
                                   cnt, datatype, op, MPI_COMM_WORLD, ierr)
                call check_ierr(ierr)
            end do
        else ! no limit
            call MPI_Allreduce(MPI_IN_PLACE, mat, size(mat), datatype, op, MPI_COMM_WORLD, ierr)
            call check_ierr(ierr)
        end if
    end subroutine allreduce_i_3

    subroutine allreduce_i_4(mat, optional_op)
        ! Allreduce for 4 dimensional integer array
        implicit none
        integer, intent(inout) :: mat(:, :, :, :)
        integer, optional, intent(in) :: optional_op
        integer :: ii, ie, ji, je, ki, ke, li, le
        integer :: i, j, k, l, idx_end, cnt
        integer :: op = op_mpi_sum ! default operation
        integer :: datatype

        if (present(optional_op)) then
            call check_operation(optional_op)
            op = optional_op
        end if

        ! Set the first and last index of the array for each dimension.
        ii = lbound(mat, 1); ie = ubound(mat, 1)
        ji = lbound(mat, 2); je = ubound(mat, 2)
        ki = lbound(mat, 3); ke = ubound(mat, 3)
        li = lbound(mat, 4); le = ubound(mat, 4)
        if (sizeof(mat(ii, ji, ki, li)) == 4) then
            datatype = MPI_INTEGER4 ! 4 byte integer
        else
            datatype = MPI_INTEGER8 ! 8 byte integer
        end if
        ! Because of the limitation of the size of the array, the array is divided into several parts and allreduce is performed.
        if (max_i4 < size(mat, 1)) then ! 1st index is larger than the limit(max_i4)
            do l = li, le
                do k = ki, ke
                    do j = ji, je
                        do i = ii, ie, max_i4
                            idx_end = min(i + max_i4 - 1, ie)
                            cnt = idx_end - i + 1
                            call MPI_Allreduce(MPI_IN_PLACE, mat(i:idx_end, j, k, l), &
                                               cnt, datatype, op, MPI_COMM_WORLD, ierr)
                            call check_ierr(ierr)
                        end do
                    end do
                end do
            end do
        else if (max_i4 < size(mat, 1)*size(mat, 2)) then ! 1st index * 2nd index is larger than the limit(max_i4)
            do l = li, le
                do k = ki, ke
                    do j = ji, je, max_i4
                        idx_end = min(j + max_i4 - 1, je)
                        cnt = (idx_end - j + 1)*size(mat, 1)
                        call MPI_Allreduce(MPI_IN_PLACE, mat(:, j:idx_end, k, l), &
                                           cnt, datatype, op, MPI_COMM_WORLD, ierr)
                        call check_ierr(ierr)
                    end do
                end do
            end do
        else if (max_i4 < size(mat, 1)*size(mat, 2)*size(mat, 3)) then ! 1st index * 2nd index * 3rd index is larger than the limit(max_i4)
            do l = li, le
                do k = ki, ke, max_i4
                    idx_end = min(k + max_i4 - 1, ke)
                    cnt = (idx_end - k + 1)*size(mat, 1)*size(mat, 2)
                    call MPI_Allreduce(MPI_IN_PLACE, mat(:, :, k:idx_end, l), &
                                       cnt, datatype, op, MPI_COMM_WORLD, ierr)
                    call check_ierr(ierr)
                end do
            end do
        else if (max_i4 < size(mat)) then ! 1st index * 2nd index * 3rd index * 4th index is larger than the limit(max_i4)
            do l = li, le, max_i4
                idx_end = min(l + max_i4 - 1, le)
                cnt = (idx_end - l + 1)*size(mat, 1)*size(mat, 2)*size(mat, 3)
                call MPI_Allreduce(MPI_IN_PLACE, mat(:, :, :, l:idx_end), &
                                   cnt, datatype, op, MPI_COMM_WORLD, ierr)
                call check_ierr(ierr)
            end do
        else ! no limit
            call MPI_Allreduce(MPI_IN_PLACE, mat, size(mat), datatype, op, MPI_COMM_WORLD, ierr)
            call check_ierr(ierr)
        end if
    end subroutine allreduce_i_4

    subroutine allreduce_r(mat, optional_op)
        ! Allreduce for a real(8) value
        implicit none
        real(8), intent(inout) :: mat
        integer, optional, intent(in) :: optional_op
        integer :: ii, ie
        integer :: i, idx_end, cnt
        integer :: op = op_mpi_sum ! default operation

        if (present(optional_op)) then
            call check_operation(optional_op)
            op = optional_op
        end if

        call MPI_Allreduce(MPI_IN_PLACE, mat, 1, MPI_REAL8, op, MPI_COMM_WORLD, ierr)
        call check_ierr(ierr)

    end subroutine allreduce_r

    subroutine allreduce_r_1(mat, optional_op)
        ! Allreduce for 1 dimensional real(8) array
        implicit none
        real(8), intent(inout) :: mat(:)
        integer, optional, intent(in) :: optional_op
        integer :: ii, ie
        integer :: i, idx_end, cnt
        integer :: op = op_mpi_sum ! default operation

        if (present(optional_op)) then
            call check_operation(optional_op)
            op = optional_op
        end if

        ! Set the first and last index of the array for each dimension.
        ii = lbound(mat, 1); ie = ubound(mat, 1)

        ! Because of the limitation of the size of the array, the array is divided into several parts and allreduce is performed.
        if (max_i4 < size(mat, 1)) then ! 1st index is larger than the limit(max_i4)
            do i = ii, ie, max_i4
                idx_end = min(i + max_i4 - 1, ie)
                cnt = idx_end - i + 1
                call MPI_Allreduce(MPI_IN_PLACE, mat(i:idx_end), &
                                   cnt, MPI_REAL8, op, MPI_COMM_WORLD, ierr)
                call check_ierr(ierr)
            end do
        else ! no limit
            call MPI_Allreduce(MPI_IN_PLACE, mat, size(mat), MPI_REAL8, op, MPI_COMM_WORLD, ierr)
            call check_ierr(ierr)
        end if
    end subroutine allreduce_r_1

    subroutine allreduce_r_2(mat, optional_op)
        ! Allreduce for 2 dimensional real(8) array
        implicit none
        real(8), intent(inout) :: mat(:, :)
        integer, optional, intent(in) :: optional_op
        integer :: ii, ie, ji, je
        integer :: i, j, idx_end, cnt
        integer :: op = op_mpi_sum ! default operation

        if (present(optional_op)) then
            call check_operation(optional_op)
            op = optional_op
        end if

        ! Set the first and last index of the array for each dimension.
        ii = lbound(mat, 1); ie = ubound(mat, 1)
        ji = lbound(mat, 2); je = ubound(mat, 2)

        ! Because of the limitation of the size of the array, the array is divided into several parts and allreduce is performed.
        if (max_i4 < size(mat, 1)) then ! 1st index is larger than the limit(max_i4)
            do j = ji, je
                do i = ii, ie, max_i4
                    idx_end = min(i + max_i4 - 1, ie)
                    cnt = idx_end - i + 1
                    call MPI_Allreduce(MPI_IN_PLACE, mat(i:idx_end, j), &
                                       cnt, MPI_REAL8, op, MPI_COMM_WORLD, ierr)
                    call check_ierr(ierr)
                end do
            end do
        else if (max_i4 < size(mat, 1)*size(mat, 2)) then ! 1st index * 2nd index is larger than the limit(max_i4)
            do j = ji, je, max_i4
                idx_end = min(j + max_i4 - 1, je)
                cnt = (idx_end - j + 1)*size(mat, 1)
                call MPI_Allreduce(MPI_IN_PLACE, mat(:, j:idx_end), &
                                   cnt, MPI_REAL8, op, MPI_COMM_WORLD, ierr)
                call check_ierr(ierr)
            end do
        else ! no limit
            call MPI_Allreduce(MPI_IN_PLACE, mat, size(mat), MPI_REAL8, op, MPI_COMM_WORLD, ierr)
            call check_ierr(ierr)
        end if
    end subroutine allreduce_r_2

    subroutine allreduce_r_3(mat, optional_op)
        ! Allreduce for 3 dimensional real(8) array
        implicit none
        real(8), intent(inout) :: mat(:, :, :)
        integer, optional, intent(in) :: optional_op
        integer :: ii, ie, ji, je, ki, ke
        integer :: i, j, k, idx_end, cnt
        integer :: op = op_mpi_sum ! default operation

        if (present(optional_op)) then
            call check_operation(optional_op)
            op = optional_op
        end if

        ! Set the first and last index of the array for each dimension.
        ii = lbound(mat, 1); ie = ubound(mat, 1)
        ji = lbound(mat, 2); je = ubound(mat, 2)
        ki = lbound(mat, 3); ke = ubound(mat, 3)

        ! Because of the limitation of the size of the array, the array is divided into several parts and allreduce is performed.
        if (max_i4 < size(mat, 1)) then ! 1st index is larger than the limit(max_i4)
            do k = ki, ke
                do j = ji, je
                    do i = ii, ie, max_i4
                        idx_end = min(i + max_i4 - 1, ie)
                        cnt = idx_end - i + 1
                        call MPI_Allreduce(MPI_IN_PLACE, mat(i:idx_end, j, k), &
                                           cnt, MPI_REAL8, op, MPI_COMM_WORLD, ierr)
                        call check_ierr(ierr)
                    end do
                end do
            end do
        else if (max_i4 < size(mat, 1)*size(mat, 2)) then ! 1st index * 2nd index is larger than the limit(max_i4)
            do k = ki, ke
                do j = ji, je, max_i4
                    idx_end = min(j + max_i4 - 1, je)
                    cnt = (idx_end - j + 1)*size(mat, 1)
                    call MPI_Allreduce(MPI_IN_PLACE, mat(:, j:idx_end, k), &
                                       cnt, MPI_REAL8, op, MPI_COMM_WORLD, ierr)
                    call check_ierr(ierr)
                end do
            end do
        else if (max_i4 < size(mat, 1)*size(mat, 2)*size(mat, 3)) then ! 1st index * 2nd index * 3rd index is larger than the limit(max_i4)
            do k = ki, ke, max_i4
                idx_end = min(k + max_i4 - 1, ke)
                cnt = (idx_end - k + 1)*size(mat, 1)*size(mat, 2)
                call MPI_Allreduce(MPI_IN_PLACE, mat(:, :, k:idx_end), &
                                   cnt, MPI_REAL8, op, MPI_COMM_WORLD, ierr)
                call check_ierr(ierr)
            end do
        else ! no limit
            call MPI_Allreduce(MPI_IN_PLACE, mat, size(mat), MPI_REAL8, op, MPI_COMM_WORLD, ierr)
            call check_ierr(ierr)
        end if
    end subroutine allreduce_r_3

    subroutine allreduce_r_4(mat, optional_op)
        ! Allreduce for 4 dimensional real(8) array
        implicit none
        real(8), intent(inout) :: mat(:, :, :, :)
        integer, optional, intent(in) :: optional_op
        integer :: ii, ie, ji, je, ki, ke, li, le
        integer :: i, j, k, l, idx_end, cnt
        integer :: op = op_mpi_sum ! default operation

        if (present(optional_op)) then
            call check_operation(optional_op)
            op = optional_op
        end if

        ! Set the first and last index of the array for each dimension.
        ii = lbound(mat, 1); ie = ubound(mat, 1)
        ji = lbound(mat, 2); je = ubound(mat, 2)
        ki = lbound(mat, 3); ke = ubound(mat, 3)
        li = lbound(mat, 4); le = ubound(mat, 4)

        ! Because of the limitation of the size of the array, the array is divided into several parts and allreduce is performed.
        if (max_i4 < size(mat, 1)) then ! 1st index is larger than the limit(max_i4)
            do l = li, le
                do k = ki, ke
                    do j = ji, je
                        do i = ii, ie, max_i4
                            idx_end = min(i + max_i4 - 1, ie)
                            cnt = idx_end - i + 1
                            call MPI_Allreduce(MPI_IN_PLACE, mat(i:idx_end, j, k, l), &
                                               cnt, MPI_REAL8, op, MPI_COMM_WORLD, ierr)
                            call check_ierr(ierr)
                        end do
                    end do
                end do
            end do
        else if (max_i4 < size(mat, 1)*size(mat, 2)) then ! 1st index * 2nd index is larger than the limit(max_i4)
            do l = li, le
                do k = ki, ke
                    do j = ji, je, max_i4
                        idx_end = min(j + max_i4 - 1, je)
                        cnt = (idx_end - j + 1)*size(mat, 1)
                        call MPI_Allreduce(MPI_IN_PLACE, mat(:, j:idx_end, k, l), &
                                           cnt, MPI_REAL8, op, MPI_COMM_WORLD, ierr)
                        call check_ierr(ierr)
                    end do
                end do
            end do
        else if (max_i4 < size(mat, 1)*size(mat, 2)*size(mat, 3)) then ! 1st index * 2nd index * 3rd index is larger than the limit(max_i4)
            do l = li, le
                do k = ki, ke, max_i4
                    idx_end = min(k + max_i4 - 1, ke)
                    cnt = (idx_end - k + 1)*size(mat, 1)*size(mat, 2)
                    call MPI_Allreduce(MPI_IN_PLACE, mat(:, :, k:idx_end, l), &
                                       cnt, MPI_REAL8, op, MPI_COMM_WORLD, ierr)
                    call check_ierr(ierr)
                end do
            end do
        else if (max_i4 < size(mat)) then ! 1st index * 2nd index * 3rd index * 4th index is larger than the limit(max_i4)
            do l = li, le, max_i4
                idx_end = min(l + max_i4 - 1, le)
                cnt = (idx_end - l + 1)*size(mat, 1)*size(mat, 2)*size(mat, 3)
                call MPI_Allreduce(MPI_IN_PLACE, mat(:, :, :, l:idx_end), &
                                   cnt, MPI_REAL8, op, MPI_COMM_WORLD, ierr)
                call check_ierr(ierr)
            end do
        else ! no limit
            call MPI_Allreduce(MPI_IN_PLACE, mat, size(mat), MPI_REAL8, op, MPI_COMM_WORLD, ierr)
            call check_ierr(ierr)
        end if
    end subroutine allreduce_r_4

    subroutine allreduce_c(mat, optional_op)
        ! Allreduce for a complex*16 value
        implicit none
        complex*16, intent(inout) :: mat
        integer, optional, intent(in) :: optional_op
        integer :: ii, ie
        integer :: i, idx_end, cnt
        integer :: op = op_mpi_sum ! default operation

        if (present(optional_op)) then
            call check_operation(optional_op)
            op = optional_op
        end if

        call MPI_Allreduce(MPI_IN_PLACE, mat, 1, MPI_COMPLEX16, op, MPI_COMM_WORLD, ierr)
        call check_ierr(ierr)

    end subroutine allreduce_c

    subroutine allreduce_c_1(mat, optional_op)
        ! Allreduce for 1 dimensional complex*16 array
        implicit none
        complex*16, intent(inout) :: mat(:)
        integer, optional, intent(in) :: optional_op
        integer :: ii, ie
        integer :: i, idx_end, cnt
        integer :: op = op_mpi_sum ! default operation

        if (present(optional_op)) then
            call check_operation(optional_op)
            op = optional_op
        end if

        ! Set the first and last index of the array for each dimension.
        ii = lbound(mat, 1); ie = ubound(mat, 1)

        ! Because of the limitation of the size of the array, the array is divided into several parts and allreduce is performed.
        if (max_i4 < size(mat, 1)) then ! 1st index is larger than the limit(max_i4)
            do i = ii, ie, max_i4
                idx_end = min(i + max_i4 - 1, ie)
                cnt = idx_end - i + 1
                call MPI_Allreduce(MPI_IN_PLACE, mat(i:idx_end), &
                                   cnt, MPI_COMPLEX16, op, MPI_COMM_WORLD, ierr)
                call check_ierr(ierr)
            end do
        else ! no limit
            call MPI_Allreduce(MPI_IN_PLACE, mat, size(mat), MPI_COMPLEX16, op, MPI_COMM_WORLD, ierr)
            call check_ierr(ierr)
        end if
    end subroutine allreduce_c_1

    subroutine allreduce_c_2(mat, optional_op)
        ! Allreduce for 2 dimensional complex*16 array
        implicit none
        complex*16, intent(inout) :: mat(:, :)
        integer, optional, intent(in) :: optional_op
        integer :: ii, ie, ji, je
        integer :: i, j, idx_end, cnt
        integer :: op = op_mpi_sum ! default operation

        if (present(optional_op)) then
            call check_operation(optional_op)
            op = optional_op
        end if

        ! Set the first and last index of the array for each dimension.
        ii = lbound(mat, 1); ie = ubound(mat, 1)
        ji = lbound(mat, 2); je = ubound(mat, 2)

        ! Because of the limitation of the size of the array, the array is divided into several parts and allreduce is performed.
        if (max_i4 < size(mat, 1)) then ! 1st index is larger than the limit(max_i4)
            do j = ji, je
                do i = ii, ie, max_i4
                    idx_end = min(i + max_i4 - 1, ie)
                    cnt = idx_end - i + 1
                    call MPI_Allreduce(MPI_IN_PLACE, mat(i:idx_end, j), &
                                       cnt, MPI_COMPLEX16, op, MPI_COMM_WORLD, ierr)
                    call check_ierr(ierr)
                end do
            end do
        else if (max_i4 < size(mat, 1)*size(mat, 2)) then ! 1st index * 2nd index is larger than the limit(max_i4)
            do j = ji, je, max_i4
                idx_end = min(j + max_i4 - 1, je)
                cnt = (idx_end - j + 1)*size(mat, 1)
                call MPI_Allreduce(MPI_IN_PLACE, mat(:, j:idx_end), &
                                   cnt, MPI_COMPLEX16, op, MPI_COMM_WORLD, ierr)
                call check_ierr(ierr)
            end do
        else ! no limit
            call MPI_Allreduce(MPI_IN_PLACE, mat, size(mat), MPI_COMPLEX16, op, MPI_COMM_WORLD, ierr)
            call check_ierr(ierr)
        end if
    end subroutine allreduce_c_2

    subroutine allreduce_c_3(mat, optional_op)
        ! Allreduce for 3 dimensional complex*16 array
        implicit none
        complex*16, intent(inout) :: mat(:, :, :)
        integer, optional, intent(in) :: optional_op
        integer :: ii, ie, ji, je, ki, ke
        integer :: i, j, k, idx_end, cnt
        integer :: op = op_mpi_sum ! default operation

        if (present(optional_op)) then
            call check_operation(optional_op)
            op = optional_op
        end if

        ! Set the first and last index of the array for each dimension.
        ii = lbound(mat, 1); ie = ubound(mat, 1)
        ji = lbound(mat, 2); je = ubound(mat, 2)
        ki = lbound(mat, 3); ke = ubound(mat, 3)

        ! Because of the limitation of the size of the array, the array is divided into several parts and allreduce is performed.
        if (max_i4 < size(mat, 1)) then ! 1st index is larger than the limit(max_i4)
            do k = ki, ke
                do j = ji, je
                    do i = ii, ie, max_i4
                        idx_end = min(i + max_i4 - 1, ie)
                        cnt = idx_end - i + 1
                        call MPI_Allreduce(MPI_IN_PLACE, mat(i:idx_end, j, k), &
                                           cnt, MPI_COMPLEX16, op, MPI_COMM_WORLD, ierr)
                        call check_ierr(ierr)
                    end do
                end do
            end do
        else if (max_i4 < size(mat, 1)*size(mat, 2)) then ! 1st index * 2nd index is larger than the limit(max_i4)
            do k = ki, ke
                do j = ji, je, max_i4
                    idx_end = min(j + max_i4 - 1, je)
                    cnt = (idx_end - j + 1)*size(mat, 1)
                    call MPI_Allreduce(MPI_IN_PLACE, mat(:, j:idx_end, k), &
                                       cnt, MPI_COMPLEX16, op, MPI_COMM_WORLD, ierr)
                    call check_ierr(ierr)
                end do
            end do
        else if (max_i4 < size(mat, 1)*size(mat, 2)*size(mat, 3)) then ! 1st index * 2nd index * 3rd index is larger than the limit(max_i4)
            do k = ki, ke, max_i4
                idx_end = min(k + max_i4 - 1, ke)
                cnt = (idx_end - k + 1)*size(mat, 1)*size(mat, 2)
                call MPI_Allreduce(MPI_IN_PLACE, mat(:, :, k:idx_end), &
                                   cnt, MPI_COMPLEX16, op, MPI_COMM_WORLD, ierr)
                call check_ierr(ierr)
            end do
        else ! no limit
            call MPI_Allreduce(MPI_IN_PLACE, mat, size(mat), MPI_COMPLEX16, op, MPI_COMM_WORLD, ierr)
            call check_ierr(ierr)
        end if
    end subroutine allreduce_c_3

    subroutine allreduce_c_4(mat, optional_op)
        ! Allreduce for 4 dimensional complex*16 array
        implicit none
        complex*16, intent(inout) :: mat(:, :, :, :)
        integer, optional, intent(in) :: optional_op
        integer :: ii, ie, ji, je, ki, ke, li, le
        integer :: i, j, k, l, idx_end, cnt
        integer :: op = op_mpi_sum ! default operation

        if (present(optional_op)) then
            call check_operation(optional_op)
            op = optional_op
        end if

        ! Set the first and last index of the array for each dimension.
        ii = lbound(mat, 1); ie = ubound(mat, 1)
        ji = lbound(mat, 2); je = ubound(mat, 2)
        ki = lbound(mat, 3); ke = ubound(mat, 3)
        li = lbound(mat, 4); le = ubound(mat, 4)

        ! Because of the limitation of the size of the array, the array is divided into several parts and allreduce is performed.
        if (max_i4 < size(mat, 1)) then ! 1st index is larger than the limit(max_i4)
            do l = li, le
                do k = ki, ke
                    do j = ji, je
                        do i = ii, ie, max_i4
                            idx_end = min(i + max_i4 - 1, ie)
                            cnt = idx_end - i + 1
                            call MPI_Allreduce(MPI_IN_PLACE, mat(i:idx_end, j, k, l), &
                                               cnt, MPI_COMPLEX16, op, MPI_COMM_WORLD, ierr)
                            call check_ierr(ierr)
                        end do
                    end do
                end do
            end do
        else if (max_i4 < size(mat, 1)*size(mat, 2)) then ! 1st index * 2nd index is larger than the limit(max_i4)
            do l = li, le
                do k = ki, ke
                    do j = ji, je, max_i4
                        idx_end = min(j + max_i4 - 1, je)
                        cnt = (idx_end - j + 1)*size(mat, 1)
                        call MPI_Allreduce(MPI_IN_PLACE, mat(:, j:idx_end, k, l), &
                                           cnt, MPI_COMPLEX16, op, MPI_COMM_WORLD, ierr)
                        call check_ierr(ierr)
                    end do
                end do
            end do
        else if (max_i4 < size(mat, 1)*size(mat, 2)*size(mat, 3)) then ! 1st index * 2nd index * 3rd index is larger than the limit(max_i4)
            do l = li, le
                do k = ki, ke, max_i4
                    idx_end = min(k + max_i4 - 1, ke)
                    cnt = (idx_end - k + 1)*size(mat, 1)*size(mat, 2)
                    call MPI_Allreduce(MPI_IN_PLACE, mat(:, :, k:idx_end, l), &
                                       cnt, MPI_COMPLEX16, op, MPI_COMM_WORLD, ierr)
                    call check_ierr(ierr)
                end do
            end do
        else if (max_i4 < size(mat)) then ! 1st index * 2nd index * 3rd index * 4th index is larger than the limit(max_i4)
            do l = li, le, max_i4
                idx_end = min(l + max_i4 - 1, le)
                cnt = (idx_end - l + 1)*size(mat, 1)*size(mat, 2)*size(mat, 3)
                call MPI_Allreduce(MPI_IN_PLACE, mat(:, :, :, l:idx_end), &
                                   cnt, MPI_COMPLEX16, op, MPI_COMM_WORLD, ierr)
                call check_ierr(ierr)
            end do
        else ! no limit
            call MPI_Allreduce(MPI_IN_PLACE, mat, size(mat), MPI_COMPLEX16, op, MPI_COMM_WORLD, ierr)
            call check_ierr(ierr)
        end if
    end subroutine allreduce_c_4

    subroutine check_operation(op)
        ! Check the operation for MPI functions
        use module_error
        implicit none
        integer, intent(in) :: op
        integer:: cnt
        cnt = count(ops == op)
        if (cnt == 0) then
            if (rank == 0) then
                print *, "MPI error: invalid operation, op_number = ", op
                print *, "MPI error: valid operations numbers are: ", ops
            end if
            call stop_with_errorcode(1)
        end if
    end subroutine check_operation

    subroutine check_ierr(err)
        use module_error
        implicit none
        integer :: err
        if (err /= 0) then
            print *, "MPI error: ierr = ", err
            call stop_with_errorcode(err)
        end if
    end subroutine check_ierr
#endif
end module module_mpi
