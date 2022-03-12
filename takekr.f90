! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE takekr(i, j, k, l, cint2)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
    integer, intent(inout)  :: i, j, k, l
    complex*16, intent(inout)  :: cint2
    real                        :: signij, signkl, signijkl

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!
! Consider Kramers pair integrals (i~j~|k~l~)*
!

    if (mod(i, 2) == 0) then; i = i - 1
    else; i = i + 1; end if

    if (mod(j, 2) == 0) then; j = j - 1
    else; j = j + 1; end if

    if (mod(k, 2) == 0) then; k = k - 1
    else; k = k + 1; end if

    if (mod(l, 2) == 0) then; l = l - 1
    else; l = l + 1; end if

    if (mod(i + j, 2) == 0) then; signij = 1.0d+00
    else; signij = -1.0d+00
    end if

    if (mod(k + l, 2) == 0) then; signkl = 1.0d+00
    else; signkl = -1.0d+00
    end if
    ! i = i - (-1)**mod(i,2)
    ! j = j - (-1)**mod(j,2)
    ! k = k - (-1)**mod(k,2)
    ! l = l - (-1)**mod(l,2)
    ! signijkl = (-1)**(mod(i+j,2) + mod(k+l,2))
    ! cint2 = signijkl*dconjg(cint2)
    cint2 = signij*signkl*DCONJG(cint2)

End subroutine takekr
