! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE solvH_ord  (e0, e2h)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE

        real*8,     intent(in) :: e0
        real*8,     intent(out):: e2h

        Integer                :: ia, ib, ii, ij, syma, sym1, sym2, i, j, k, l
        Integer                :: i0, j0, tab, nab, tij, nij, count
        Integer,allocatable    :: ia0(:), ib0(:), ii0(:), ij0(:), iab(:,:), iij(:,:)
        Complex*16             :: cint2
        Complex*16,allocatable :: v(:,:)
        Real*8                 :: e, signij, signkl
        Integer                :: iii, jjj

        real*8  :: thresd

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE H IS NOW CALCULATED
!      
!     EaiEbj|0>         a > b, i > j
!
!   DRAS1 =-2   DRAS2 = 0    DRAS3 = +2 
!
!
!  S(ckdl,aibj) = d(ac)d(bd)d(lj)d(ik)
!
!  (H0-E0S)ckdl,aibj = d(ac)d(bd)d(lj)d(ik)(eps(a)+eps(b)-eps(i)-eps(j)) = e(a,b,i,j)
!
!  V(aibj)   = (ai|bj) - (aj|bi)
!   
! E2h = V(aibj)/e(a,b,i,j)


!        thresd = 1.0D-08
!        thres = 1.0D-08

        e2h = 0.0d+00
        e = 0.0d+00

        i0 = 0
        Do ia = ninact+nact+1, ninact+nact+nsec
           Do ib = ninact+nact+1, ia-1
              i0 = i0 + 1
           Enddo
        Enddo

        nab = i0

        Allocate(iab(ninact+nact+1:ninact+nact+nsec,ninact+nact+1:ninact+nact+nsec))
        Allocate(ia0(nab))
        Allocate(ib0(nab))
        iab = 0

        i0 = 0
        Do ia = ninact+nact+1, ninact+nact+nsec
           Do ib = ninact+nact+1, ia-1
              i0 = i0 + 1
              iab(ia, ib) = i0
              iab(ib, ia) = i0
              ia0(i0) = ia
              ib0(i0) = ib
           Enddo
        Enddo

        i0 = 0
        Do ii =1,  ninact
           Do ij =1,  ii-1
              i0 = i0 + 1
           Enddo
        Enddo

        nij = i0
        Allocate(iij(1:ninact,1:ninact))
        Allocate(ii0(nij))
        Allocate(ij0(nij))
        iij = 0
        
        i0 = 0
        Do ii =1,  ninact
           Do ij =1,  ii-1
              i0 = i0 + 1
              iij(ii, ij) = i0
              iij(ij, ii) = i0
              ii0(i0) = ii
              ij0(i0) = ij
           Enddo
        Enddo

        Allocate (v(nab, nij))
        v = 0.0d+00

        open(1, file ='Hint', status='old', form='unformatted')
 30     read(1, err=10, end=20) i,j,k,l,cint2
        count = 0

 40     if(i<=k .or. j==l) goto 30

!         write(*,*)i,j,k,l,cint2

        tab = iab(i, k)
        tij = iij(j, l)

!        write(*,*)tab,iab(i,k),i,k
!  V(aibj)   = (ai|bj) - (aj|bi)     i > j, a > b

        if    ( i > k .and. j > l) then
           v(tab, tij) = v(tab, tij) + cint2
        
        elseif( i > k .and. j < l) then
           v(tab, tij) = v(tab, tij) - cint2

        elseif( i < k .and. j > l) then               ! (kl|ij)  l > j + ; l < j -
           v(tab, tij) = v(tab, tij) - cint2

        elseif( i < k .and. j < l) then
           v(tab, tij) = v(tab, tij) + cint2

        endif

        goto 30

!Iwamuro modify
        Do iii = 1, tab
          Do jjj = 1, tij
             IF(abs(v(iii, jjj)) > 1.0E-08)then
               write(*,'("i,j,V_h  ",2I4,2E15.7)') iii, jjj, v(iii, jjj)
             Endif
          Enddo
        Enddo

 20     close(1)

        write(*,*)'reading int2 is over'

        Do i0 = 1, nab
           ia = ia0(i0)
           ib = ib0(i0)
           sym1 = MULTB(irpmo(ia), nsymrpa+1)
           sym1 = MULTB(irpmo(ib), sym1)
           Do j0 = 1, nij
              ii = ii0(j0)
              ij = ij0(j0)
              sym2 = MULTB2(irpmo(ii), sym1)
              sym2 = MULTB2(irpmo(ij), sym2)
!Iwamuro modify
!            sym1 = MULTB2(irpmo(ia), nsymrpa+1)
!            sym1 = MULTB(irpmo(ii), sym1)
!            sym2 = MULTB2(irpmo(ij),nsymrpa+1)
!            sym2 = MULTB(irpmo(ib), sym2)

              if(sym2 == nsymrpa+1) then
!Iwamuro modify 
!              if(sym1 == sym2) then

                 e = eps(ia) + eps(ib) - eps(ii) - eps(ij) + eshift  ! For Level Shift (2007/2/9)

                 coeff1 = v(i0, j0)/e
                 sumc2local = sumc2local + ABS(coeff1)**2

                 e2h = e2h - DCONJG(v(i0, j0))*v(i0, j0)/e
              endif
           End do
        End do

        write(*,'("e2h      = ",E20.10,"a.u.")')e2h

        write(*,'("sumc2,h  = ",E20.10)')sumc2local
        sumc2 = sumc2 + sumc2local



        deallocate(v)
        deallocate(iab)
        deallocate(ia0)
        deallocate(ib0)
        deallocate(iij)
        deallocate(ii0)
        deallocate(ij0)


 10     continue                !write(*,*)'error about opening Hint file' ;stop
 100    continue 
      write(*,*)'end solvh_ord'
   End SUBROUTINE solvH_ord





      