! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE tramo1 ( i, j, int1)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE
        integer,    intent(in)  :: i, j
        complex*16, intent(out) :: int1


        integer :: i0, j0, sym1, sym2
        integer :: n(2,2), mo(2)


! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=


     int1  = 0.0d+00
     n(:,:)= 0
     sym1 = irpamo(i)
     sym2 = irpamo(j)

     If(sym1 == sym2) then

        mo(1) = i
        mo(2) = j

        Do i0 = 1, 2
           if( mo(i0) <= ninact ) then
              n(i0,1) = 1 
              n(i0,2) = ninact
           elseif( mo(i0) >= ninact+1 .and. mo(i0) <= ninact+nact ) then
              n(i0,1) = ninact+1 
              n(i0,2) = ninact+nact
           elseif( mo(i0) >= ninact+nact+1 .and. mo(i0) <= ninact+nact+nsec ) then
              n(i0,1) = ninact+nact+1 
              n(i0,2) = ninact+nact+nsec
           endif
        End do    ! i0

        do i0 = n(1,1), n(1,2)
        do j0 = n(2,1), n(2,2)
           If(irpamo(i0) ==sym1 .and. irpamo(j0) ==sym2) then
              int1 = int1 + DCONJG(f(i0,i))*CMPLX(oner(i0,j0),onei(i0,j0),16)*f(j0,j)
           Endif
        end do
        end do

      Endif

        End subroutine tramo1




! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE tramo2 ( i, j, k, l, int2)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE
        integer,    intent(in) :: i, j, k, l

        complex*16,intent(out) :: int2


        integer     :: i0, j0, k0, l0, sym1, sym2, sym3, sym4, sym5, sym6
        integer     :: n(4,2), mo(4)
        integer     :: nint, tcount, count
                    
        real*8      :: i2r, i2i, nsign
        complex*16  :: cmplxint

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=


     int2     = 0.0d+00
     cmplxint = 0.0d+00
     n    = 0
     sym1 = irpamo(i)
     sym2 = irpamo(j)
     sym3 = irpamo(k)
     sym4 = irpamo(l)

     sym5 = MULTB2(sym1,nsymrpa+1)
     sym5 = MULTB (sym2, sym5)
     sym6 = MULTB2(sym4,nsymrpa+1) 
     sym6 = MULTB (sym3, sym6)

!           if (i==1 .and. j==1.and. k==2.and.l==2) then
!             write(*,*)"sym",sym1,sym2,sym3,sym4,sym5,sym6 
!           endif

!     If(MULTB(sym1,sym2) == MULTB(sym3,sym4)) then
     If(sym5 == sym6) then

!Iwamuro modify
!           if (i==1 .and. j==1.and. k==2.and.l==2) then
!             write(*, *) "1122" 
!           endif

        mo(1) = i
        mo(2) = j
        mo(3) = k
        mo(4) = l

        Do i0 = 1, 4

           if( mo(i0) <= ninact ) then

              n(i0,1) = 1 
              n(i0,2) = ninact

           elseif( mo(i0) >= ninact+1 .and. mo(i0) <= ninact+nact ) then

              n(i0,1) = ninact+1 
              n(i0,2) = ninact+nact

           elseif( mo(i0) >= ninact+nact+1 .and. mo(i0) <= ninact+nact+nsec ) then

              n(i0,1) = ninact+nact+1
              n(i0,2) = ninact+nact+nsec


           endif
! Iwamuro modify
!       write(*,*) " mo(i0), n(i0,1), n(i0,2) =", mo(i0), n(i0,1), n(i0,2)
!        if(debug) write(*,*) mo(i0), n(i0, 1), n(i0, 2)

        End do    ! i0
        tcount = 0
        count  = 0


        do i0 = n(1,1), n(1,2)
        do j0 = n(2,1), n(2,2)
        do k0 = n(3,1), n(3,2)
        do l0 = n(4,1), n(4,2)
           tcount = tcount + 1

         If(irpamo(i0) == sym1 .and. irpamo(j0) == sym2 .and. &
         &  irpamo(k0) == sym3 .and. irpamo(l0) == sym4 ) then

!           if (i==1 .and. j==1.and. k==2.and.l==2) then
!             write(*, '("1122",4I4,10E15.5)') i0, j0, k0, l0, i2r, i2i, (f(i0,i)),(f(k0,k)),(j0,j),f(l0,l)
!           endif

!           if (i==1 .and. j==1.and. k==1.and.l==1) then
!             write(*, '("1111",4I4,10E15.5)') i0, j0, k0, l0, i2r, i2i, (f(i0,i)),(f(k0,k)),f(j0,j),f(l0,l)
!           endif  

         count = count +1
           cmplxint = 0.0d+00

           nint = ABS(indtwr(i0,j0,k0,l0))
           nsign = SIGN(1,indtwr(i0,j0,k0,l0))
            i2r = int2r(nint)*nsign
           nsign = SIGN(1,indtwi(i0,j0,k0,l0))
           i2i = int2i(nint)*nsign

           cmplxint = CMPLX(i2r, i2i, 16)
           
            int2 = int2 + DCONJG(f(i0,i))*DCONJG(f(k0,k))*f(j0,j)*f(l0,l)*cmplxint

!Iwamuro modify
!         write(*,*) "nint, nsign, i2r =", nint, nsign, i2r
!         write(*,*) "i2i, cmplxint =", i2i, cmplxint
!         write(*,'(4I4, 2E15.5)') i0, j0, k0, l0, i2r, i2i
!         write(*,'(4I4,8E15.5)') i0, j0, k0, l0, f(i0,i),f(j0,j),f(k0,k),f(l0,l)
!           write(*,'(4I4,4E15.5)') i0, j0, k0, l0, int2, cmplxint
!         write(*,'(8E10.4)') f(i0,i),f(j0,j),f(k0,k),f(l0,l)
         


         Endif

        end do
        end do
        end do
        end do

! Iwamuro modify
!        write(*,*) "tcount, count =", tcount, count

      Endif

      End subroutine tramo2

