module module_takekr

    implicit none

    private
    public :: takekr

    interface takekr
        module procedure takekr_complex, takekr_real
    end interface takekr
contains
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE takekr_complex(i, j, k, l, cint2)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: sign_even_ret1, sign_odd_ret1

        Implicit NONE
        integer, intent(inout)      :: i, j, k, l
        complex*16, intent(inout)   :: cint2
        real                        :: signij, signkl

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!
! Consider Kramers pair integrals (i~j~|k~l~)*
!

        i = i + sign_odd_ret1(i) ! i = i+1 if i is odd, otherwise i = i-1
        j = j + sign_odd_ret1(j)
        k = k + sign_odd_ret1(k)
        l = l + sign_odd_ret1(l)

        signij = sign_even_ret1(i + j) ! signij = 1 if i+j is even, otherwise signij = -1
        signkl = sign_even_ret1(k + l)

        cint2 = signij*signkl*DCONJG(cint2)

    End subroutine takekr_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE takekr_real(i, j, k, l, rint2)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: sign_even_ret1, sign_odd_ret1

        Implicit NONE
        integer, intent(inout)      :: i, j, k, l
        real(8), intent(inout)      :: rint2
        real                        :: signij, signkl

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!
! Consider Kramers pair integrals (i~j~|k~l~)*
!

        i = i + sign_odd_ret1(i) ! i = i+1 if i is odd, otherwise i = i-1
        j = j + sign_odd_ret1(j)
        k = k + sign_odd_ret1(k)
        l = l + sign_odd_ret1(l)

        signij = sign_even_ret1(i + j) ! signij = 1 if i+j is even, otherwise signij = -1
        signkl = sign_even_ret1(k + l)

        rint2 = signij*signkl*rint2

    End subroutine takekr_real

end module module_takekr
