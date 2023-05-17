! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE dets(fa, occold, occnew, ds)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables
    use module_index_utils, only: sign_even_ret1

    Implicit NONE
    real(8), intent(in)  :: fa(ninact + 1:ninact + nact, ninact + 1:ninact + nact)
    integer, intent(in)  :: occold(nelec), occnew(nelec)
    complex*16, intent(out) :: ds

    integer                  :: i0, j0, i, info, ipvt(nelec), n
    complex*16               :: sini(nelec, nelec), work(nelec), det(2), ds2
    real(8)                   :: phase

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!       Calculating determinant of S(k|k~)

    sini = 0.0d+00
    work = 0.0d+00
    det = 0.0d+00
    ds = 0.0d+00

    Do i0 = 1, nelec
        Do j0 = 1, nelec
            sini(i0, j0) = fa(occold(i0) + ninact, occnew(j0) + ninact)
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

End subroutine dets

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE detsc(fac, occold, occnew, ds)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables
    use module_index_utils, only: sign_even_ret1

    Implicit NONE
    complex*16, intent(in)  :: fac(ninact + 1:ninact + nact, ninact + 1:ninact + nact)
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
            sini(i0, j0) = fac(occold(i0) + ninact, occnew(j0) + ninact)
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
