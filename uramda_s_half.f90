! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE uramda_s_half (ur, wnew, dimn, dimm) 

! Assume C1 molecule, calculate U(ramda_S)^(-1/2) : ramda_S is diagonal matrix of S
!                                                   but eliminated 0 eigen values
!                                                   the dimension beocmes dimm


! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE
     
        integer, intent(in)    :: dimn, dimm           ! dimn >= dimm reduced dimension
        real*8, intent(inout)  :: ur(dimn,dimm)
        real*8, intent(in)     :: wnew(dimm)

        real*8              :: ramdah(dimm,dimm)                  ! ramdaS^(-1/2)  dimm*dimm matrix
        integer             :: j0, j, i, i0, i1
        integer             :: k0, l0, ii, jj, kk, ll


        ramdah(:,:) = 0.0d+00

        do i0 = 1, dimm
           ramdah(i0,i0) = wnew(i0)**(-5.0d-01)
        end do

        ur = MATMUL(ur,ramdah) ! ur: N*M ramdah: M*M  then (ur)(randah) : N*M

!        write(*,*) 'urams'
!        do i0 = 1, dimn
!        do j0 = 1, dimm
!           if (ABS(urams(i0,j0))>=0.0001) then
!              write(*,*) urams (i0,j0)
!           endif
!        end do
!        end do

   end subroutine uramda_s_half



! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE ucramda_s_half (uc, wnew, dimn, dimm) 

! Assume C1 molecule, calculate U(ramda_S)^(-1/2) : ramda_S is diagonal matrix of S
!                                                   but eliminated 0 eigen values
!                                                   the dimension beocmes dimm


! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE
     
        integer, intent(in)        :: dimn, dimm           ! dimn >= dimm reduced dimension
        complex*16, intent(inout)  :: uc(dimn,dimm)
        real*8, intent(in)         :: wnew(dimm)

        real*8              :: ramdah(dimm,dimm)                  ! ramdaS^(-1/2)  dimm*dimm matrix
        integer             :: j0, j, i, i0, i1
        integer             :: k0, l0, ii, jj, kk, ll


        ramdah(:,:) = 0.0d+00

        do i0 = 1, dimm
           ramdah(i0,i0) = wnew(i0)**(-5.0d-01)
        end do

        uc = MATMUL(uc,ramdah) ! ur: N*M ramdah: M*M  then (ur)(randah) : N*M

!        write(*,*) 'uc'
!        do i0 = 1, dimn
!        do j0 = 1, dimm
!           if (ABS(uc(i0,j0))>=0.0001) then
!              write(*,'("uc",2E15.8,2I4)') uc(i0,j0),i0,j0
!           endif
!        end do
!        end do

   end subroutine ucramda_s_half



