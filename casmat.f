! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   SUBROUTINE casmat(mat)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        use four_caspt2_module

        Implicit NONE

        complex*16, intent(out) :: mat(ndet, ndet)


        integer              :: occ, vir, indr, inds, inda, indb
        integer              :: ir, is, ia, ib, imo, nint
        integer              :: i0, j0, k0, l0, i, j, newidet1, newidet2
        integer              :: phase, phase1, phase2
        real*8               :: nsign, i2r, i2i
        complex*16           :: cmplxint, mat0

        integer, allocatable :: ridet(:), oc(:), vi(:)

        mat = 0.0d+00

        write(*,*)'Cas mat enter'

        Allocate(oc(nelec))
        Allocate(vi(nact-nelec))

        Do i = 1, ndet

           occ = 0
           oc  = 0
           vir = 0
           vi  = 0

           Do imo = 1, nact
              If(BTEST(idet(i), imo-1)) then
                 occ = occ + 1
                 oc(occ) = imo
              Else 
                 vir = vir + 1
                 vi(vir) = imo
              Endif
           Enddo

!           write(*,*) 'i, idet(i)',i, idet(i)
!           write(*,*) occ, oc(1:occ)
!           write(*,*) vir, vi(1:vir)
!           write(*,*) ' '


!! IDENTICAL DETERMINANT => DIAGONAL TERM 

           cmplxint = 0.0d+00

           Do i0 = 1, ninact
              ir = i0
              cmplxint = CMPLX(oner(ir,ir), onei(ir,ir), 16)
              mat(i, i) = mat(i, i) + cmplxint
           End do

           Do i0 = 1, nelec
              indr = oc(i0)
              ir = indr+ninact
              cmplxint = CMPLX(oner(ir,ir), onei(ir,ir), 16)
              mat(i, i) = mat(i, i) + cmplxint
           End do

           mat0 = 0.0d+00

           Do i0 = 1, ninact+nelec

              if(i0 <= ninact) then
                 ir = i0
              Else
                 indr = i0 - ninact
                 indr = oc(indr)
                 ir   = indr + ninact
              End if
              Do j0 = i0+1, ninact+nelec

                 if(j0 <= ninact) then
                    is = j0
                 Else
                    inds = j0 - ninact
                    inds = oc(inds)
                    is   = inds + ninact
                 End if


                 nint = ABS(indtwr(ir,ir,is,is))
                 nsign = SIGN(1,indtwr(ir,ir,is,is))
                 i2r = int2r(nint)*nsign
                 nsign = SIGN(1,indtwi(ir,ir,is,is))
                 i2i = int2i(nint)*nsign
                 cmplxint = CMPLX(i2r, i2i, 16)
!                 write(*,*)ir,is,cmplxint

                 mat0 = mat0 + 0.5d+00*cmplxint
!                 mat(i,i) = mat(i,i) + 0.5d+00*cmplxint

                 nint = ABS(indtwr(ir,is,is,ir))
                 nsign = SIGN(1,indtwr(ir,is,is,ir))
                 i2r = int2r(nint)*nsign
                 nsign = SIGN(1,indtwi(ir,is,is,ir))
                 i2i = int2i(nint)*nsign
                 cmplxint = CMPLX(i2r, i2i, 16)
!                 write(*,*)ir,is,cmplxint

                 mat0 = mat0 - 0.5d+00*cmplxint
!                 mat(i,i) = mat(i,i) - 0.5d+00*cmplxint

              End do
           End do

           mat(i,i) = mat(i,i) + mat0 + DCONJG(mat0)

!           write(*,*)'mat(',i,',',i,') = ',mat(i,i)

!! ONE SPINOR DIFFERENCE

           Do i0 = 1, nelec
              indr = oc(i0)
              ir = indr+ninact

              Do k0 = 1, nact - nelec
                 inda = vi(k0)
                 ia = inda+ninact

                 Call one_e_exct (idet(i), inda, indr, newidet1, phase1)

                 Call idetr(newidet1, j)

!                 write(*,*)'j=',j

                 If(j > i) then


                    cmplxint = CMPLX(oner(ir,ia), onei(ir,ia), 16)
                    mat(i, j) = mat(i, j) + cmplxint


                    Do l0 = 1, ninact
                       is = l0

                       nint = ABS(indtwr(ir,ia,is,is))
                       nsign = SIGN(1,indtwr(ir,ia,is,is))
                       i2r = int2r(nint)*nsign
                       nsign = SIGN(1,indtwi(ir,ia,is,is))
                       i2i = int2i(nint)*nsign
                       cmplxint = CMPLX(i2r, i2i, 16)

                       mat(i, j) = mat(i, j) + cmplxint
                 
                       nint = ABS(indtwr(ir,is,is,ia))
                       nsign = SIGN(1,indtwr(ir,is,is,ia))
                       i2r = int2r(nint)*nsign
                       nsign = SIGN(1,indtwi(ir,is,is,ia))
                       i2i = int2i(nint)*nsign
                       cmplxint = CMPLX(i2r, i2i, 16)

                       mat(i, j) = mat(i, j) - cmplxint
                    End do      !l0


                    Do l0 = 1, nelec
                       inds = oc(l0)
                       is = inds+ninact

                       nint = ABS(indtwr(ir,ia,is,is))
                       nsign = SIGN(1,indtwr(ir,ia,is,is))
                       i2r = int2r(nint)*nsign
                       nsign = SIGN(1,indtwi(ir,ia,is,is))
                       i2i = int2i(nint)*nsign
                       cmplxint = CMPLX(i2r, i2i, 16)

                       mat(i, j) = mat(i, j) + cmplxint
                 
                       nint = ABS(indtwr(ir,is,is,ia))
                       nsign = SIGN(1,indtwr(ir,is,is,ia))
                       i2r = int2r(nint)*nsign
                       nsign = SIGN(1,indtwi(ir,is,is,ia))
                       i2i = int2i(nint)*nsign
                       cmplxint = CMPLX(i2r, i2i, 16)

                       mat(i, j) = mat(i, j) - cmplxint
                    End do      !l0

                    if( mod(phase1, 2) == 0) phase =  1.0d+00
                    if( mod(phase1, 2) == 1) phase = -1.0d+00
                    
                    mat(i, j) = phase * mat(i, j)
                    mat(j, i) = DCONJG(mat(i, j))

                 Endif

              End do            ! k0
           End do               ! i0

!! TWO ELECTRON DIFFERNT CASE

           Do i0 = 1, nelec
              Do j0 = i0 + 1, nelec
                 indr = oc(i0)
                 inds = oc(j0)
                 ir = indr+ninact
                 is = inds+ninact
!                 write(*,*)'ir,indr,is,inds',ir,indr,is,inds

                 Do k0 = 1, nact - nelec
                    Do l0 = k0 + 1, nact - nelec
                       inda = vi(k0)
                       indb = vi(l0)
                       ia = inda+ninact
                       ib = indb+ninact


                       Call one_e_exct (idet(i), inda, indr, newidet1, phase1)
                       Call one_e_exct (newidet1, indb, inds, newidet2, phase2)

                       Call idetr(newidet2, j)

                       If(j > i) then
                       
                          if( mod(phase1 + phase2, 2) == 0) phase =  1.0d+00
                          if( mod(phase1 + phase2, 2) == 1) phase = -1.0d+00

                          nint = ABS(indtwr(ir,ia,is,ib))
                          nsign = SIGN(1,indtwr(ir,ia,is,ib))
                          i2r = int2r(nint)*nsign
                          nsign = SIGN(1,indtwi(ir,ia,is,ib))
                          i2i = int2i(nint)*nsign
                          cmplxint = CMPLX(i2r, i2i, 16)

                          mat(i, j) = cmplxint
                 
                          nint = ABS(indtwr(ir,ib,is,ia))
                          nsign = SIGN(1,indtwr(ir,ib,is,ia))
                          i2r = int2r(nint)*nsign
                          nsign = SIGN(1,indtwi(ir,ib,is,ia))
                          i2i = int2i(nint)*nsign
                          cmplxint = CMPLX(i2r, i2i, 16)

                          mat(i, j) = mat(i, j) - cmplxint

                          mat(i, j) = phase * mat(i, j)
                          mat(j, i) = DCONJG(mat(i, j))

                       Endif

                    End do      ! l0
                 End do         ! k0
              End do            ! j0
           End do               ! i0

        End do                  ! i

        Deallocate(oc)
        Deallocate(vi)

 1000   end subroutine


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   SUBROUTINE idetr(iidet, j)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        use four_caspt2_module

        Implicit NONE

        integer, intent(in)  :: iidet
        integer, intent(out) :: j
        integer              :: i0

        j = 0

        Do i0 = 1, ndet
           If(idet(i0) == iidet) then
              j = i0
           Endif
        End do

   END SUBROUTINE idetr
