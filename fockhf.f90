! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE fockhf ! TO CALCULATE FOCK MATRIX OF HF STATE, A TEST

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
    integer :: ii, jj, kk, ll
    integer :: j, i, k, l
    integer :: nint, n

    real*8 :: i2r, i2i, dr, di, nsign
    complex*16 :: cmplxint, dens

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    debug = .TRUE.
    thres = 1.0d-15
!        thres = 0.0d+00

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

!! TEST TO CALCULATE FOCK MATRIX OF HF STATE fpq = hpq + SIGUMA_r[(pq|rr)-(pr|qr)]
!! THIS MUST BE DIAGONAL MATRIX AND DIAGONAL ELEMENTS CORESPONDS TO SPINOR ENERGIES.
    write (*, *) ' '
    write (*, *) 'FOR TEST, FOCK MATRIX OF HF STATE IS CALCULATED '

    n = 0
    f = 0.0d+00

    do i = 1, ninact + nact + nsec
        do j = i, ninact + nact + nsec

            f(i, j) = DCMPLX(oner(i, j), onei(i, j))

            do k = 1, ninact + nelec

                i2r = 0.0d+00
                i2i = 0.0d+00
                dr = 0.0d+00
                di = 0.0d+00
                cmplxint = 0.0d+00

                nint = ABS(indtwr(i, j, k, k))

                nsign = SIGN(1, indtwr(i, j, k, k))
                i2r = int2r(nint)*nsign

                nsign = SIGN(1, indtwi(i, j, k, k))
                i2i = int2i(nint)*nsign

                cmplxint = CMPLX(i2r, i2i, 16)

                nint = ABS(indtwr(i, k, k, j))

                nsign = SIGN(1, indtwr(i, k, k, j))
                i2r = int2r(nint)*nsign

                nsign = SIGN(1, indtwi(i, k, k, j))
                i2i = int2i(nint)*nsign

                cmplxint = cmplxint - CMPLX(i2r, i2i, 16)

                f(i, j) = f(i, j) + cmplxint
!                     write(*,*)f(i,j)
            End do           ! k

            f(j, i) = DCONJG(f(i, j))

        End do       ! j
    End do          ! i

    write (*, *) ' '
    write (*, *) 'OFF DIAGONAL ELEMENTS OF FOCK MATRIX WHICH IS LARGER THAN 1.0d-06 '
    write (*, *) ' '
    do i = 1, ninact + nact + nsec
        do j = i, ninact + nact + nsec
            if ((i /= j) .and. (ABS(f(i, j)) > 1.0d-6)) then
!            if(i/=j)then
                write (*, '(2I4,2E20.10)') i, j, f(i, j)
            end if
        end do
    end do
    write (*, *) ' '
    write (*, *) 'THESE DIAGONAL ELEMENTS SHOULD BE CORESPOND TO HF SPINOR ENERGY '
    write (*, *) ' '
    write (*, *) '  NO.   Spinor Energy(Re)   Spinor Energy(Im) '&
    &, 'Spinor Energy (HF)        ERROR'
    do i = 1, ninact + nact + nsec
        write (*, '(I4,4E20.10)') i, f(i, i), orbmo(i), orbmo(i) - dble(f(i, i))
    end do

    write (*, *) 'fockhf end'

end SUBROUTINE fockhf
