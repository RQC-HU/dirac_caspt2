! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE tramo1_ty(i, j, int1)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE
       integer, intent(in)  :: i, j
       complex*16, intent(out) :: int1

       integer :: i0, j0, sym1, sym2
       integer :: n(2, 2), mo(2)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       int1 = 0.0d+00
       n(:, :) = 0
       sym1 = irpamo(i)
       sym2 = irpamo(j)

       If (MULTB_D(sym1, sym2) == 1) then

           mo(1) = i
           mo(2) = j

           Do i0 = 1, 2
               if (mo(i0) <= ninact) then
                   n(i0, 1) = 1
                   n(i0, 2) = ninact
               elseif (mo(i0) >= ninact + 1 .and. mo(i0) <= ninact + nact) then
                   n(i0, 1) = ninact + 1
                   n(i0, 2) = ninact + nact
               elseif (mo(i0) >= ninact + nact + 1 .and. mo(i0) <= ninact + nact + nsec) then
                   n(i0, 1) = ninact + nact + 1
                   n(i0, 2) = ninact + nact + nsec
               end if
           End do    ! i0

           do i0 = n(1, 1), n(1, 2)
           do j0 = n(2, 1), n(2, 2)
               If (irpamo(i0) == sym1 .and. irpamo(j0) == sym2) then
                  !! Adding one-electron integral to the fock matrics is executed only by the master process
                  !! because DIRAC's one-electron integral file (MRCONEE) is not
                  !! devided even if DIRAC is executed in parallel (MPI).
                   if (rank == 0) then
                       int1 = int1 + DCONJG(f(i0, i))*CMPLX(oner(i0, j0), onei(i0, j0), 16)*f(j0, j)
                   end if
               End if
           end do
           end do

       End if

   End subroutine tramo1_ty

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE tramo2_ty(i, j, k, l, int2)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE
       integer, intent(in) :: i, j, k, l

       complex*16, intent(out) :: int2

       integer     :: i0, j0, k0, l0, sym1, sym2, sym3, sym4, sym5, sym6
       integer     :: n(4, 2), mo(4)
       integer     :: nint, tcount, count

       real*8      :: i2r, i2i, nsign, nsign2
       complex*16  :: cmplxint

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       int2 = 0.0d+00
       cmplxint = 0.0d+00
       n = 0
       sym1 = irpamo(i)
       sym2 = irpamo(j)
       sym3 = irpamo(k)
       sym4 = irpamo(l)

       sym5 = MULTB_D(sym1, sym2)
       sym6 = MULTB_D(sym3, sym4)

       ! Iwamuro modify
!       write(*,*) "sym1 =", sym1
!       write(*,*) "sym2 =", sym2
!       write(*,*) "sym3 =", sym3
!       write(*,*) "sym4 =", sym4
!       write(*,*) "sym5 =", sym5
!       write(*,*) "sym6 =", sym6

       If (MULTB_S(sym5, sym6) == 1) then

           mo(1) = i
           mo(2) = j
           mo(3) = k
           mo(4) = l

           Do i0 = 1, 4

               if (mo(i0) <= ninact) then

                   n(i0, 1) = 1
                   n(i0, 2) = ninact

               elseif (mo(i0) >= ninact + 1 .and. mo(i0) <= ninact + nact) then

                   n(i0, 1) = ninact + 1
                   n(i0, 2) = ninact + nact

               elseif (mo(i0) >= ninact + nact + 1 .and. mo(i0) <= ninact + nact + nsec) then

                   n(i0, 1) = ninact + nact + 1
                   n(i0, 2) = ninact + nact + nsec

               end if

! Iwamuro modify
!       write(*,*) " mo(i0), n(i0,1), n(i0,2) =", mo(i0), n(i0,1), n(i0,2)
!        if(debug) write(*,*) mo(i0), n(i0, 1), n(i0, 2)

           End do    ! i0
           tcount = 0
           count = 0

           do i0 = n(1, 1), n(1, 2)
           do j0 = n(2, 1), n(2, 2)
           do k0 = n(3, 1), n(3, 2)
           do l0 = n(4, 1), n(4, 2)
               tcount = tcount + 1

               If (irpamo(i0) == sym1 .and. irpamo(j0) == sym2 .and. &
               &  irpamo(k0) == sym3 .and. irpamo(l0) == sym4) then
                   count = count + 1
                   cmplxint = 0.0d+00

                   nint = ABS(indtwr(i0, j0, k0, l0))
                   nsign = SIGN(1, indtwr(i0, j0, k0, l0))
                   i2r = int2r(nint)*nsign
                   nsign2 = SIGN(1, indtwi(i0, j0, k0, l0))
                   i2i = int2i(nint)*nsign2

                   cmplxint = CMPLX(i2r, i2i, 16)

                   int2 = int2 + DCONJG(f(i0, i))*DCONJG(f(k0, k))*f(j0, j)*f(l0, l)*cmplxint

!Iwamuro modify
!         write(*,*) "nint, nsign, i2r =", nint, nsign, i2r
!         write(*,*) "nsign2, i2i, cmplxint =", nsign2, i2i, cmplxint
!         write(*,*) "int2 =", int2

               End if

           end do
           end do
           end do
           end do

!Iwamuro modify
!        write(*,*) "tcount, count =", tcount, count
!        write(*,*)tcount, count
       End if

   End subroutine tramo2_ty
