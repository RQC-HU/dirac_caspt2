module module_ulambda_s_half

    implicit none
    private
    public ulambda_s_half
    interface ulambda_s_half
        module procedure ulambda_s_half_real, ulambda_s_half_complex
    end interface ulambda_s_half

contains
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE ulambda_s_half_real(ur, wnew, dimn, dimm)

! Assume C1 molecule, calculate U(ramda_S)^(-1/2) : ramda_S is diagonal matrix of S
!                                                   but eliminated 0 eigen values
!                                                   the dimension beocmes dimm

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use four_caspt2_module

        Implicit NONE

        integer, intent(in)    :: dimn, dimm           ! dimn >= dimm reduced dimension
        real*8, intent(inout)  :: ur(dimn, dimm)
        real*8, intent(in)     :: wnew(dimm)

        real*8              :: ramdah(dimm, dimm)                  ! ramdaS^(-1/2)  dimm*dimm matrix
        integer             :: i0

        ramdah(:, :) = 0.0d+00

        do i0 = 1, dimm
            ramdah(i0, i0) = wnew(i0)**(-5.0d-01)
        end do

        ur = MATMUL(ur, ramdah) ! ur: N*M ramdah: M*M  then (ur)(randah) : N*M

    end subroutine ulambda_s_half_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE ulambda_s_half_complex(uc, wnew, dimn, dimm)

! Assume C1 molecule, calculate U(ramda_S)^(-1/2) : ramda_S is diagonal matrix of S
!                                                   but eliminated 0 eigen values
!                                                   the dimension beocmes dimm

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use four_caspt2_module

        Implicit NONE

        integer, intent(in)        :: dimn, dimm           ! dimn >= dimm reduced dimension
        complex*16, intent(inout)  :: uc(dimn, dimm)
        real*8, intent(in)         :: wnew(dimm)

        real*8              :: ramdah(dimm, dimm)                  ! ramdaS^(-1/2)  dimm*dimm matrix
        integer             :: i0

        ramdah(:, :) = 0.0d+00

        do i0 = 1, dimm
            ramdah(i0, i0) = wnew(i0)**(-5.0d-01)
        end do

        uc = MATMUL(uc, ramdah) ! ur: N*M ramdah: M*M  then (ur)(randah) : N*M

    end subroutine ulambda_s_half_complex
end module module_ulambda_s_half
