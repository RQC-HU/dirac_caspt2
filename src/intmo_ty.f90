
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE intmo1_ty(i, j, int1)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
    integer, intent(in)  :: i, j
    complex*16, intent(out) :: int1

    integer :: sym1, sym2

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    sym1 = irpamo(i)
    sym2 = irpamo(j)

    If (MULTB_D(sym1, sym2) == 1) then
        int1 = CMPLX(oner(i, j), onei(i, j), 16)
!     else
!       int1 = 0.0d+00
    End if

End subroutine intmo1_ty

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE intmo2_ty(i, j, k, l, int2)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
    integer, intent(in) :: i, j, k, l

    complex*16, intent(out) :: int2

    integer     :: sym1, sym2, sym3, sym4, syma, symb, symc
    integer     :: nint

    real*8      :: i2r, i2i, nsign

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    sym1 = irpamo(i)
    sym2 = irpamo(j)
    sym3 = irpamo(k)
    sym4 = irpamo(l)
    syma = MULTB_D(sym1, sym2)
    symb = MULTB_D(sym3, sym4)
    symc = MULTB_S(syma, symb)

    If (symc == 1) then

        !   nint = ABS(indtwr(i, j, k, l))
        !   nsign = SIGN(1,indtwr(i, j, k, l))
        !   i2r = int2r(nint)*nsign
        i2r = inttwr(i, j, k, l)
        !   nsign = SIGN(1,indtwi(i, j, k, l))
        !   i2i = int2i(nint)*nsign
        i2i = inttwi(i, j, k, l)

        int2 = CMPLX(i2r, i2i, 16)

!Iwamuro modify
!        write(*,*)'nint', nint
!        write(*,*)'nsign', nsign
!        write(*,*)'i2r', i2r
!        write(*,'(3E15.5)')int2
!
    else
        int2 = 0.0d+00
    End if

End subroutine intmo2_ty
