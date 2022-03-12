! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE intmo1(i, j, int1)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE
       integer, intent(in)  :: i, j
       complex*16, intent(out) :: int1

       integer :: sym1, sym2

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       int1 = 0.0d+00
       sym1 = irpamo(i)
       sym2 = irpamo(j)

!     If(sym1 == sym2) then
       int1 = CMPLX(oner(i, j), onei(i, j), 16)
!     End if

   End subroutine intmo1

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE intmo2(i, j, k, l, int2)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE
       integer, intent(in) :: i, j, k, l

       complex*16, intent(out) :: int2

       integer     :: sym1, sym2, sym3, sym4
       integer     :: nint

       real*8      :: i2r, i2i, nsign

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       int2 = 0.0d+00
       sym1 = irpamo(i)
       sym2 = irpamo(j)
       sym3 = irpamo(k)
       sym4 = irpamo(l)

       !    If(MULTB(sym1,sym2) == MULTB(sym3,sym4)) then

       if (i == 15 .and. j == 3 .and. k == 4 .and. l == 4) write (*, *) 'int number', ABS(indtwr(i, j, k, l))

       nint = ABS(indtwr(i, j, k, l))
       nsign = SIGN(1, indtwr(i, j, k, l))
       i2r = int2r(nint)*nsign
       nsign = SIGN(1, indtwi(i, j, k, l))
       i2i = int2i(nint)*nsign

       int2 = CMPLX(i2r, i2i, 16)

!Iwamuro modify
!        write(*,'(3E15.5)')int2

!      Endif

   End subroutine intmo2
