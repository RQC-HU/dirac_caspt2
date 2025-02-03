module module_blas
    implicit none

    interface gemv
        module procedure dgemv_
        module procedure zgemv_
    end interface gemv

    interface gemm
        module procedure dgemm_
        module procedure zgemm_
    end interface gemm
contains

    ! vec_ret = mat_a * vec_b (real, 2D matrix x 1D vector)
    subroutine dgemv_(mat_a, vec_b, vec_ret)
        implicit none
        real*8, intent(in) :: mat_a(:, :)
        real*8, intent(in) :: vec_b(:)
        real*8, intent(out) :: vec_ret(size(mat_a, 1))
        integer :: m, n, lda

        vec_ret = 0.0d0
        m = size(mat_a, 1)
        n = size(mat_a, 2)
        lda = size(mat_a, 1)
        call dgemv('N', m, n, 1.0d0, mat_a, lda, vec_b, 1, 0.0d0, vec_ret, 1)

    end subroutine dgemv_

    ! vec_ret = mat_a * vec_b (complex, 2D matrix x 1D vector)
    subroutine zgemv_(mat_a, vec_b, vec_ret)
        implicit none
        complex*16, intent(in) :: mat_a(:, :)
        complex*16, intent(in) :: vec_b(:)
        complex*16, intent(out) :: vec_ret(size(mat_a, 1))
        integer :: m, n, lda

        vec_ret = 0.0d0
        m = size(mat_a, 1)
        n = size(mat_a, 2)
        lda = size(mat_a, 1)
        call zgemv('N', m, n, (1.0d0, 0.0d0), mat_a, lda, vec_b, 1, (0.0d0, 0.0d0), vec_ret, 1)
    end subroutine zgemv_

    ! mat_ret = mat_a * mat_b (real, 2D x 2D matrix)
    subroutine dgemm_(mat_a, mat_b, mat_ret)
        implicit none
        real*8, intent(in) :: mat_a(:, :)
        real*8, intent(in) :: mat_b(:, :)
        real*8, intent(out) :: mat_ret(size(mat_a, 1), size(mat_b, 2))
        integer :: m, n, k, lda, ldb, ldc

        mat_ret = 0.0d+00
        m = size(mat_a, 1)
        n = size(mat_b, 2)
        k = size(mat_a, 2)
        lda = size(mat_a, 1)
        ldb = size(mat_b, 1)
        ldc = size(mat_a, 1)
        call dgemm('N', 'N', m, n, k, 1.0d0, mat_a, lda, mat_b, ldb, 0.0d0, mat_ret, ldc)
    end subroutine dgemm_

    ! mat_ret = mat_a * mat_b (complex, 2D x 2D matrix)
    subroutine zgemm_(mat_a, mat_b, mat_ret)
        implicit none
        complex*16, intent(in) :: mat_a(:, :)
        complex*16, intent(in) :: mat_b(:, :)
        complex*16, intent(out) :: mat_ret(size(mat_a, 1), size(mat_b, 2))
        integer :: m, n, k, lda, ldb, ldc

        mat_ret = 0.0d+00
        m = size(mat_a, 1)
        n = size(mat_b, 2)
        k = size(mat_a, 2)
        lda = size(mat_a, 1)
        ldb = size(mat_b, 1)
        ldc = size(mat_a, 1)
        call zgemm('N', 'N', m, n, k, (1.0d0, 0.0d0), mat_a, lda, mat_b, ldb, (0.0d0, 0.0d0), mat_ret, ldc)
    end subroutine zgemm_

end module module_blas
