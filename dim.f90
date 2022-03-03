! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine dim1_density(creat1, anhi1, sr, si)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       use four_caspt2_module

       Implicit NONE
       integer, intent(in) :: creat1, anhi1
       real*8, intent(out) :: sr, si

       integer :: j0, j, i, i0, i1
       integer :: k0, l0, ii, jj, kk, ll
       complex*16 :: cmplxcii, cmplxcij, cmplxs

       integer ::  newidet, phase, phasenew, nbitsa

! calculation of <0|Ec1a1|0>

       sr = 0.0d+00
       si = 0.0d+00

       do i0 = 1, ndet

           i = idet(i0)

           call one_e_exct(i, creat1, anhi1, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phase
           j0 = 0

           do i1 = 1, ndet
               j = idet(i1)
               if (j == i) then
                   j0 = i1
                   goto 1
               end if
           end do
1          continue

           if (j0 == 0) then
               go to 10
           end if

!        write(*,*)creat1,anhi1,i0,j0,phase
!        write(*,*)cir(i0,iroot) , cir(j0,iroot) , cii(i0,iroot) ,cii(j0,iroot)

!        write(*,*) 'i0,j0,phase ,cir(i0,iroot) , cir(j0,iroot) , cii(i0,iroot) ,cii(j0,iroot)'

           cmplxcii = CMPLX(cir(i0, iroot), cii(i0, iroot), 16)
           cmplxcij = CMPLX(cir(j0, iroot), cii(j0, iroot), 16)

!        cmplxs = cmplxcij*DCONJG(cmplxcii)
           cmplxs = cmplxcii*DCONJG(cmplxcij)

           if (mod(phasenew, 2) == 0) then
               sr = sr + REAL(cmplxs, 8)
               si = si + DIMAG(cmplxs)
           else
               sr = sr - REAL(cmplxs, 8)
               si = si - DIMAG(cmplxs)
           end if

!        if(creat1==1.and.anhi1==1) write(*,*)1,1,sr,si

!        if(mod(phasenew,2)==0) then
!           sr = sr + cir(i0,iroot) * cir(j0,iroot) + cii(i0,iroot) * cii(j0,iroot)
!           si = si + cii(i0,iroot) * cir(j0,iroot) - cir(i0,iroot) * cii(j0,iroot)
!        else
!           sr = sr - cir(i0,iroot) * cir(j0,iroot) + cii(i0,iroot) * cii(j0,iroot)
!           si = si - cii(i0,iroot) * cir(j0,iroot) - cir(i0,iroot) * cii(j0,iroot)
!        end if
!
!        write(*,*)'sr',sr
!        write(*,*)'si',si
10     end do

   end subroutine dim1_density

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine dim1_density_nondiag(creat1, anhi1, s)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       use four_caspt2_module

       Implicit NONE
       integer, intent(in)     :: creat1, anhi1
       complex*16, intent(out) :: s

       integer :: j0, j, i, i0, i1
       integer :: k0, l0, ii, jj, kk, ll
       complex*16 :: cmplxcii, cmplxcij, cmplxs

       integer ::  newidet, phase, phasenew, nbitsa

       s = 0.0d+00
       do i0 = 1, ndet
           i = idet(i0)
           call one_e_exct(i, creat1, anhi1, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phase
           j0 = 0

           do i1 = 1, ndet
               j = idet(i1)
               if (j == i) then
                   j0 = i1
                   goto 1
               end if
           end do
1          continue

           if (j0 == 0) then
               go to 10
           end if

           cmplxcii = CMPLX(cir(i0, iroot), cii(i0, iroot), 16)
           cmplxcij = CMPLX(cir(j0, iroot), cii(j0, iroot), 16)

           cmplxs = DCONJG(cmplxcii)*DCONJG(cmplxcij)    ! This part is unique for Aperp term.

           if (mod(phasenew, 2) == 0) then
               s = s + cmplxs
           else
               s = s - cmplxs
           end if

10     end do

   end subroutine dim1_density_nondiag

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine dim1_density_diag(creat1, anhi1, s)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       use four_caspt2_module

       Implicit NONE
       integer, intent(in)     :: creat1, anhi1
       complex*16, intent(out) :: s

       integer :: j0, j, i, i0, i1
       integer :: k0, l0, ii, jj, kk, ll
       complex*16 :: cmplxcii, cmplxcij, cmplxs

       integer ::  newidet, phase, phasenew, nbitsa

       s = 0.0d+00
       do i0 = 1, ndet
           i = idet(i0)
           call one_e_exct(i, creat1, anhi1, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phase
           j0 = 0

           do i1 = 1, ndet
               j = idet(i1)
               if (j == i) then
                   j0 = i1
                   goto 1
               end if
           end do
1          continue

           if (j0 == 0) then
               go to 10
           end if

           cmplxcii = CMPLX(cir(i0, iroot), cii(i0, iroot), 16)
           cmplxcij = CMPLX(cir(j0, iroot), cii(j0, iroot), 16)

           cmplxs = cmplxcii*DCONJG(cmplxcij)

           if (mod(phasenew, 2) == 0) then
               s = s + cmplxs
           else
               s = s - cmplxs
           end if

10     end do

   end subroutine dim1_density_diag

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine dim2_density(creat1, anhi1, creat2, anhi2, sr, si)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       use four_caspt2_module

       Implicit NONE
       integer, intent(in) :: creat2, anhi2, anhi1, creat1
       real*8, intent(out) :: sr, si
       integer :: newidet, phase, phasenew, nbitsa
       integer :: j0, j, i, i0, i1
       integer :: k0, l0, ii, jj, kk, ll
       complex*16 :: cmplxcii, cmplxcij, cmplxs

! calculation of <0|Ec1a1Ec2a2|0>

       sr = 0.0d+00
       si = 0.0d+00

       do i0 = 1, ndet

           i = idet(i0)

           call one_e_exct(i, creat2, anhi2, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phase

           call one_e_exct(i, creat1, anhi1, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phasenew + phase

           j0 = 0
           do i1 = 1, ndet
               j = idet(i1)
               if (j == i) then
                   j0 = i1
                   goto 1
               end if
           end do
1          continue

           if (j0 == 0) then
               go to 10
           end if

!        write(*,*)creat1,anhi1,i0,j0,phase
!        write(*,*)cir(i0,iroot) , cir(j0,iroot) , cii(i0,iroot) ,cii(j0,iroot)

!       Caluculation of C(i,iroot)*conjugate(C(j,iroot))

!        write(*,*) 'i0,j0,phasenew ,cir(i0,iroot) , cir(j0,iroot) , cii(i0,iroot) ,cii(j0,iroot)'
!        write(*,*) i0,j0,phasenew ,cir(i0,iroot) , cir(j0,iroot) , cii(i0,iroot) ,cii(j0,iroot)

           cmplxcii = CMPLX(cir(i0, iroot), cii(i0, iroot), 16)
           cmplxcij = CMPLX(cir(j0, iroot), cii(j0, iroot), 16)

!        cmplxs = cmplxcij*DCONJG(cmplxcii)
           cmplxs = cmplxcii*DCONJG(cmplxcij)

           if (mod(phasenew, 2) == 0) then
               sr = sr + DBLE(cmplxs)
               si = si + DIMAG(cmplxs)
           else
               sr = sr - DBLE(cmplxs)
               si = si - DIMAG(cmplxs)
           end if

!        if(mod(phasenew,2)==0) then
!           sr = sr + cir(i0,iroot) * cir(j0,iroot) + cii(i0,iroot) * cii(j0,iroot)
!           si = si + cii(i0,iroot) * cir(j0,iroot) - cir(i0,iroot) * cii(j0,iroot)
!        else
!           sr = sr - cir(i0,iroot) * cir(j0,iroot) + cii(i0,iroot) * cii(j0,iroot)
!           si = si - cii(i0,iroot) * cir(j0,iroot) - cir(i0,iroot) * cii(j0,iroot)
!        end if
!
!!        sr = sr + (-1) **phasenew *(cir(i0,iroot)*cir(j0,iroot) + cii(i0,iroot)*cii(j0,iroot))
!!                 ! (cir(i0,iroot)*cir(j0,iroot)+(-1)(- cii(j0,iroot))*cii(i0,iroot))
!!
!!        si = si + (-1) **phasenew *(cii(i0,iroot)*cir(j0,iroot) - cir(i0,iroot)*cii(j0,iroot))
!!                 ! (cii(i0,iroot)*cir(j0,iroot)+(- cii(j0,iroot))*cir(i0,iroot))
!!!        write(*,*)'sr',sr
!!!        write(*,*)'si',si

10     end do

   end subroutine dim2_density

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine dim3_density(creat1, anhi1, creat2, anhi2, creat3, anhi3, sr, si)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       use four_caspt2_module

       Implicit NONE
       integer, intent(in) :: anhi3, anhi2, anhi1, creat3, creat2, creat1
       real*8, intent(out) :: sr, si
       integer :: newidet, phase, phasenew, nbitsa
       integer :: j0, j, i, i0, i1
       integer :: k0, l0, ii, jj, kk, ll
       complex*16 :: cmplxcii, cmplxcij, cmplxs

! calculation of <0|Ec1a1Ec2a2Ec3a3|0>

       sr = 0.0d+00
       si = 0.0d+00

       do i0 = 1, ndet

           i = idet(i0)

           call one_e_exct(i, creat3, anhi3, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phase

           call one_e_exct(i, creat2, anhi2, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phasenew + phase

           call one_e_exct(i, creat1, anhi1, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phasenew + phase

           j0 = 0
           do i1 = 1, ndet
               j = idet(i1)
               if (j == i) then
                   j0 = i1
                   goto 1
               end if
           end do
1          continue

           if (j0 == 0) then
               go to 10
           end if

!        write(*,*) 'i0,j0,phasenew ,cir(i0,iroot) , cir(j0,iroot) , cii(i0,iroot) ,cii(j0,iroot)'
!        write(*,*) i0,j0,phasenew ,cir(i0,iroot) , cir(j0,iroot) , cii(i0,iroot) ,cii(j0,iroot)

           cmplxcii = CMPLX(cir(i0, iroot), cii(i0, iroot), 16)
           cmplxcij = CMPLX(cir(j0, iroot), cii(j0, iroot), 16)

!        cmplxs = cmplxcij*DCONJG(cmplxcii)
           cmplxs = cmplxcii*DCONJG(cmplxcij)

           if (mod(phasenew, 2) == 0) then
               sr = sr + DBLE(cmplxs)
               si = si + DIMAG(cmplxs)
           else
               sr = sr - DBLE(cmplxs)
               si = si - DIMAG(cmplxs)
           end if

!        if(mod(phasenew,2)==0) then
!           sr = sr + cir(i0,iroot) * cir(j0,iroot) + cii(i0,iroot) * cii(j0,iroot)
!           si = si + cii(i0,iroot) * cir(j0,iroot) - cir(i0,iroot) * cii(j0,iroot)
!        else
!           sr = sr - cir(i0,iroot) * cir(j0,iroot) + cii(i0,iroot) * cii(j0,iroot)
!           si = si - cii(i0,iroot) * cir(j0,iroot) - cir(i0,iroot) * cii(j0,iroot)
!        end if

!!        sr = sr + (-1) **phasenew *(cir(i0,iroot)*cir(j0,iroot) + cii(i0,iroot)*cii(j0,iroot))
!!                 ! (cir(i0,iroot)*cir(j0,iroot)+(-1)(- cii(j0,iroot))*cii(i0,iroot))
!!
!!        si = si + (-1) **phasenew *(cii(i0,iroot)*cir(j0,iroot) - cir(i0,iroot)*cii(j0,iroot))
!!                 ! (cii(i0,iroot)*cir(j0,iroot)+(- cii(j0,iroot))*cii(i0,iroot))

10     end do

   end subroutine dim3_density

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine dim4_density(creat1, anhi1, creat2, anhi2, creat3, anhi3, creat4, anhi4, sr, si)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       use four_caspt2_module

       Implicit NONE
       integer, intent(in) :: anhi4, anhi3, anhi2, anhi1, creat4, creat3, creat2, creat1
       real*8, intent(out) :: sr, si
       integer :: newidet, phase, phasenew, nbitsa
       integer :: j0, j, i, i0, i1
       integer :: k0, l0, ii, jj, kk, ll
       complex*16 :: cmplxcii, cmplxcij, cmplxs

! calculation of <0|Ec1a1Ec2a2Ec3a3Ec4a4|0>

       sr = 0.0d+00
       si = 0.0d+00

       do i0 = 1, ndet

           i = idet(i0)

           call one_e_exct(i, creat4, anhi4, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phase

           call one_e_exct(i, creat3, anhi3, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phasenew + phase

           call one_e_exct(i, creat2, anhi2, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phasenew + phase

           call one_e_exct(i, creat1, anhi1, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet

           phasenew = phasenew + phase
           j0 = 0
           do i1 = 1, ndet
               j = idet(i1)
               if (j == i) then
                   j0 = i1
                   goto 1
               end if
           end do
1          continue

           if (j0 == 0) then
               go to 10
           end if

!        write(*,*) 'i0,j0,phasenew ,cir(i0,iroot) , cir(j0,iroot) , cii(i0,iroot) ,cii(j0,iroot)'
!        write(*,*) i0,j0,phasenew ,cir(i0,iroot) , cir(j0,iroot) , cii(i0,iroot) ,cii(j0,iroot)

           cmplxcii = CMPLX(cir(i0, iroot), cii(i0, iroot), 16)
           cmplxcij = CMPLX(cir(j0, iroot), cii(j0, iroot), 16)

!        cmplxs = cmplxcij*DCONJG(cmplxcii)
           cmplxs = cmplxcii*DCONJG(cmplxcij)

           if (mod(phasenew, 2) == 0) then
               sr = sr + DBLE(cmplxs)
               si = si + DIMAG(cmplxs)
           else
               sr = sr - DBLE(cmplxs)
               si = si - DIMAG(cmplxs)
           end if

!        if(mod(phasenew,2)==0) then
!           sr = sr + cir(i0,iroot) * cir(j0,iroot) + cii(i0,iroot) * cii(j0,iroot)
!           si = si + cii(i0,iroot) * cir(j0,iroot) - cir(i0,iroot) * cii(j0,iroot)
!        else
!           sr = sr - cir(i0,iroot) * cir(j0,iroot) + cii(i0,iroot) * cii(j0,iroot)
!           si = si - cii(i0,iroot) * cir(j0,iroot) - cir(i0,iroot) * cii(j0,iroot)
!        end if

!!        sr = sr + (-1) **phasenew *(cir(i0,iroot)*cir(j0,iroot) + cii(i0,iroot)*cii(j0,iroot))
!!                 ! (cir(i0,iroot)*cir(j0,iroot)+(-1)(- cii(j0,iroot))*cii(i0,iroot))
!!
!!        si = si + (-1) **phasenew *(cii(i0,iroot)*cir(j0,iroot) - cir(i0,iroot)*cii(j0,iroot))
!!                 ! (cii(i0,iroot)*cir(j0,iroot)+(- cii(j0,iroot))*cii(i0,iroot))

10     end do

   end subroutine dim4_density

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine dim1_density_R(creat1, anhi1, sr)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       use four_caspt2_module

       Implicit NONE
       integer, intent(in) :: creat1, anhi1
       real*8, intent(out) :: sr
       integer :: j0, j, i, i0, i1
       integer :: k0, l0, ii, jj, kk, ll

       integer ::  newidet, phase, phasenew, nbitsa

! calculation of <0|Ec1a1|0> for REAL

       sr = 0.0d+00

       do i0 = 1, ndet

           i = idet(i0)

           call one_e_exct(i, creat1, anhi1, newidet, phase)
           if (newidet == 0) goto 10
!        write(*,*) i, newidet, phase, creat1-ninact, anhi1-ninact

           i = newidet
           phasenew = phase
           j0 = 0

           do i1 = 1, ndet
               j = idet(i1)
               if (j == i) then
                   j0 = i1
                   goto 1
               end if
           end do
1          continue

           if (j0 == 0) then
               go to 10
           end if

!        write(*,*) 'i0,j0,iroot,phase ,cir(i0,ir)cir(j0,ir)cii(i0,ir)cii(j0,ir)'
!        write(*,*) i0,j0,iroot,phase &
!        , cir(i0,iroot) , cir(j0,iroot) , cii(i0,iroot) ,cii(j0,iroot)

           if (mod(phasenew, 2) == 0) then
               sr = sr + cir(i0, iroot)*cir(j0, iroot)
           else
               sr = sr - cir(i0, iroot)*cir(j0, iroot)
           end if

!         sr = sr + (-1)**phasenew * cir(i0,iroot) * cir(j0,iroot)

10     end do

   end subroutine dim1_density_R

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine dim2_density_R(creat1, anhi1, creat2, anhi2, sr)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       use four_caspt2_module

       Implicit NONE
       integer, intent(in) :: creat2, anhi2, anhi1, creat1
       real*8, intent(out) :: sr
       integer :: newidet, phase, phasenew, nbitsa
       integer :: j0, j, i, i0, i1
       integer :: k0, l0, ii, jj, kk, ll

! calculation of <0|Ec1a1Ec2a2|0>

       sr = 0.0d+00

       do i0 = 1, ndet

           i = idet(i0)

           call one_e_exct(i, creat2, anhi2, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phase

           call one_e_exct(i, creat1, anhi1, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phasenew + phase

           j0 = 0
           do i1 = 1, ndet
               j = idet(i1)
               if (j == i) then
                   j0 = i1
                   goto 1
               end if
           end do
1          continue

           if (j0 == 0) then
               go to 10
           end if

!       Caluculation of C(i,iroot)*conjugate(C(j,iroot))

!        write(*,*) 'i0,j0,phasenew ,cir(i0,iroot) , cir(j0,iroot) , cii(i0,iroot) ,cii(j0,iroot)'
!        write(*,*) i0,j0,phasenew ,cir(i0,iroot) , cir(j0,iroot) , cii(i0,iroot) ,cii(j0,iroot)

           if (mod(phasenew, 2) == 0) then
               sr = sr + cir(i0, iroot)*cir(j0, iroot)
           else
               sr = sr - cir(i0, iroot)*cir(j0, iroot)
           end if

!!        sr = sr + (-1) **phasenew *cir(i0,iroot)*cir(j0,iroot)

10     end do

   end subroutine dim2_density_R

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine dim3_density_R(creat1, anhi1, creat2, anhi2, creat3, anhi3, sr)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       use four_caspt2_module

       Implicit NONE
       integer, intent(in) :: anhi3, anhi2, anhi1, creat3, creat2, creat1
       real*8, intent(out) :: sr
       integer :: newidet, phase, phasenew, nbitsa
       integer :: j0, j, i, i0, i1
       integer :: k0, l0, ii, jj, kk, ll

! calculation of <0|Ec1a1Ec2a2Ec3a3|0>

       sr = 0.0d+00

       do i0 = 1, ndet

           i = idet(i0)

           call one_e_exct(i, creat3, anhi3, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phase

           call one_e_exct(i, creat2, anhi2, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phasenew + phase

           call one_e_exct(i, creat1, anhi1, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phasenew + phase

           j0 = 0
           do i1 = 1, ndet
               j = idet(i1)
               if (j == i) then
                   j0 = i1
                   goto 1
               end if
           end do
1          continue

           if (j0 == 0) then
               go to 10
           end if

!        write(*,*) 'i0,j0,phasenew ,cir(i0,iroot) , cir(j0,iroot) , cii(i0,iroot) ,cii(j0,iroot)'
!        write(*,*) i0,j0,phasenew ,cir(i0,iroot) , cir(j0,iroot) , cii(i0,iroot) ,cii(j0,iroot)

           if (mod(phasenew, 2) == 0) then
               sr = sr + cir(i0, iroot)*cir(j0, iroot)
           else
               sr = sr - cir(i0, iroot)*cir(j0, iroot)
           end if

!!         sr = sr + (-1) **phasenew * cir(i0,iroot)*cir(j0,iroot)

10     end do

   end subroutine dim3_density_R

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine dim4_density_R(creat1, anhi1, creat2, anhi2, creat3, anhi3, creat4, anhi4, sr)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       use four_caspt2_module

       Implicit NONE
       integer, intent(in) :: anhi4, anhi3, anhi2, anhi1, creat4, creat3, creat2, creat1
       real*8, intent(out) :: sr
       integer :: newidet, phase, phasenew, nbitsa
       integer :: j0, j, i, i0, i1
       integer :: k0, l0, ii, jj, kk, ll

! calculation of <0|Ec1a1Ec2a2Ec3a3Ec4a4|0>

       sr = 0.0d+00

       do i0 = 1, ndet

           i = idet(i0)

           call one_e_exct(i, creat4, anhi4, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phase

           call one_e_exct(i, creat3, anhi3, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phasenew + phase

           call one_e_exct(i, creat2, anhi2, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet
           phasenew = phasenew + phase

           call one_e_exct(i, creat1, anhi1, newidet, phase)
           if (newidet == 0) goto 10
           i = newidet

           phasenew = phasenew + phase
           j0 = 0
           do i1 = 1, ndet
               j = idet(i1)
               if (j == i) then
                   j0 = i1
                   goto 1
               end if
           end do
1          continue

           if (j0 == 0) then
               go to 10
           end if

!        write(*,*) 'i0,j0,phasenew ,cir(i0,iroot) , cir(j0,iroot) , cii(i0,iroot) ,cii(j0,iroot)'
!        write(*,*) i0,j0,phasenew ,cir(i0,iroot) , cir(j0,iroot) , cii(i0,iroot) ,cii(j0,iroot)

           if (mod(phasenew, 2) == 0) then
               sr = sr + cir(i0, iroot)*cir(j0, iroot)
           else
               sr = sr - cir(i0, iroot)*cir(j0, iroot)
           end if

!!         sr = sr + (-1) **phasenew *cir(i0,iroot)*cir(j0,iroot)

10     end do

   end subroutine dim4_density_R
