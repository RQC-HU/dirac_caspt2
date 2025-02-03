module module_takekr

    iuse, intrinsic :: iso_fortran_env, only: int32, int64
    implicit none

    private
    public :: takekr

    interface takekr
        module procedure takekr_int32_complex, takekr_int64_complex, &
                         takekr_int32_real, takekr_int64_real
    end interface takekr
contains
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE takekr_int32_complex(i, j, k, l, cint2)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables

        Implicit NONE
        integer(kind=int32), intent(inout)      :: i, j, k, l
        complex*16, intent(inout)   :: cint2
        integer                     :: signij, signkl

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!
! Consider Kramers pair integrals (i~j~|k~l~)*
!

        i = i + merge(1, -1, mod(i, 2)==1) ! i = i+1 if i is odd, otherwise i = i-1
        j = j + merge(1, -1, mod(j, 2)==1)
        k = k + merge(1, -1, mod(k, 2)==1)
        l = l + merge(1, -1, mod(l, 2)==1)

        signij = merge(1, -1, mod(i + j, 2)==0) ! signij = 1 if i+j is even, otherwise signij = -1
        signkl = merge(1, -1, mod(k + l, 2)==0)

        cint2 = signij*signkl*DCONJG(cint2)

    End subroutine takekr_int32_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE takekr_int64_complex(i, j, k, l, cint2)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables

        Implicit NONE
        integer(kind=int64), intent(inout)      :: i, j, k, l
        complex*16, intent(inout)   :: cint2
        integer                     :: signij, signkl

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!
! Consider Kramers pair integrals (i~j~|k~l~)*
!

        i = i + merge(1, -1, mod(i, 2)==1) ! i = i+1 if i is odd, otherwise i = i-1
        j = j + merge(1, -1, mod(j, 2)==1)
        k = k + merge(1, -1, mod(k, 2)==1)
        l = l + merge(1, -1, mod(l, 2)==1)

        signij = merge(1, -1, mod(i + j, 2)==0) ! signij = 1 if i+j is even, otherwise signij = -1
        signkl = merge(1, -1, mod(k + l, 2)==0)

        cint2 = signij*signkl*DCONJG(cint2)

    End subroutine takekr_int64_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE takekr_int32_real(i, j, k, l, rint2)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables

        Implicit NONE
        integer(kind=int32), intent(inout)      :: i, j, k, l
        real(8), intent(inout)      :: rint2
        integer                     :: signij, signkl

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!
! Consider Kramers pair integrals (i~j~|k~l~)*
!

        i = i + merge(1, -1, mod(i, 2)==1) ! i = i+1 if i is odd, otherwise i = i-1
        j = j + merge(1, -1, mod(j, 2)==1)
        k = k + merge(1, -1, mod(k, 2)==1)
        l = l + merge(1, -1, mod(l, 2)==1)

        signij = merge(1, -1, mod(i + j, 2)==0) ! signij = 1 if i+j is even, otherwise signij = -1
        signkl = merge(1, -1, mod(k + l, 2)==0)

        rint2 = signij*signkl*rint2

    End subroutine takekr_int32_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE takekr_int64_real(i, j, k, l, rint2)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables

        Implicit NONE
        integer(kind=int64), intent(inout)      :: i, j, k, l
        real(8), intent(inout)      :: rint2
        integer                     :: signij, signkl

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!
! Consider Kramers pair integrals (i~j~|k~l~)*
!

        i = i + merge(1, -1, mod(i, 2)==1) ! i = i+1 if i is odd, otherwise i = i-1
        j = j + merge(1, -1, mod(j, 2)==1)
        k = k + merge(1, -1, mod(k, 2)==1)
        l = l + merge(1, -1, mod(l, 2)==1)

        signij = merge(1, -1, mod(i + j, 2)==0) ! signij = 1 if i+j is even, otherwise signij = -1
        signkl = merge(1, -1, mod(k + l, 2)==0)

        rint2 = signij*signkl*rint2

    End subroutine takekr_int64_real

end module module_takekr
