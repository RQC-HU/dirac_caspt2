! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE dets(fa, occold, occnew, ds)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables
    use module_index_utils, only: convert_active_to_global_idx, sign_even_ret1

    Implicit NONE
    real(8), intent(in)  :: fa(global_act_start:global_act_end, global_act_start:global_act_end)
    integer, intent(in)  :: occold(nelec), occnew(nelec)
    real(8), intent(out) :: ds

    integer :: i0, j0, i, info, ipvt(nelec), n
    real(8) :: sini(nelec, nelec), work(nelec), det(2)
    real(8) :: phase

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!       Calculating determinant of S(k|k~)

    sini = 0.0d+00
    work = 0.0d+00
    det = 0.0d+00

    Do i0 = 1, nelec
        Do j0 = 1, nelec
            sini(i0, j0) = fa(convert_active_to_global_idx(occold(i0)), convert_active_to_global_idx(occnew(j0)))
        End do
    End do

    ! LU factorization by lapack routine dgetrf
    ! dgetrf: https://netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga0019443faea08275ca60a734d0593e60.html
    ! sini is the matrix to be factorized
    ! sini is overwritten by the LU factorization of sini
    ! ipvt is the pivot indices of the LU factorization.
    Call dgetrf(nelec, nelec, sini, nelec, ipvt, info)   ! From lapack LU fatorization!

    ! count the number of ipvt(i) /= i
    n = count(ipvt /= [(i, i=1, nelec)])
    phase = sign_even_ret1(n) ! If n is even, phase = 1.0d+00, else phase = -1.0d+00

    ds = 1.0d+00
    Do i = 1, nelec
        ds = sini(i, i)*ds
    End do
    ds = ds*phase

End subroutine dets

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE detsc(fac, occold, occnew, ds)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables
    use module_index_utils, only: convert_active_to_global_idx , sign_even_ret1

    Implicit NONE
    complex*16, intent(in)  :: fac(global_act_start:global_act_end, global_act_start:global_act_end)
    integer, intent(in)  :: occold(nelec), occnew(nelec)
    complex*16, intent(out) :: ds

    integer                  :: i0, j0, i, info, ipvt(nelec), n
    complex*16               :: sini(nelec, nelec), work(nelec), det(2), ds2
    real(8)                   :: phase

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!       Calculating determinant of S(k|k~)

!       Set elements

    sini = 0.0d+00
    work = 0.0d+00
    det = 0.0d+00
    ds = 0.0d+00

    Do i0 = 1, nelec
        Do j0 = 1, nelec
            sini(i0, j0) = fac(convert_active_to_global_idx(occold(i0)), convert_active_to_global_idx(occnew(j0)))
        End do
    End do

    Call zgetrf(nelec, nelec, sini, nelec, ipvt, info)   ! From lapack LU fatorization!

    ! count the number of ipvt(i) /= i
    n = count(ipvt /= [(i, i=1, nelec)])

    phase = sign_even_ret1(n) ! If n is even, phase = 1.0d+00, else phase = -1.0d+00

    ds2 = 1.0d+00

    Do i = 1, nelec
        ds2 = sini(i, i)*ds2
    End do

    ds2 = ds2*phase

    ds = ds2

End subroutine detsc
