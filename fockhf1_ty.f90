! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE fockhf1_ty ! TO CALCULATE FOCK MATRIX OF HF STATE, A TEST

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
       if (rank == 0) then ! Process limits for output
           write (3000, *) ' '
           write (3000, *) 'FOR TEST, FOCK MATRIX OF HF STATE IS CALCULATED '
       end if
       n = 0
       f = 0.0d+00

       do i = 1, ninact + nact
           do j = i, ninact + nact
            !! Adding one-electron integral to the fock matrics is executed only by the master process
            !! because DIRAC's one-electron integral file (MRCONEE) is not
            !! devided even if DIRAC is executed in parallel (MPI).
               if (rank == 0) then
                   f(i, j) = DCMPLX(oner(i, j), onei(i, j))
               end if
               do k = 1, ninact + nelec

                   Call intmo2_ty(i, j, k, k, cmplxint)

                   f(i, j) = f(i, j) + cmplxint

                   Call intmo2_ty(i, k, k, j, cmplxint)

                   f(i, j) = f(i, j) - cmplxint

!iwamuro modify
!                     write(*,*)f(i,j)
               End do           ! k

               f(j, i) = DCONJG(f(i, j))
           End do       ! j
       End do          ! i

       do i = ninact + nact + 1, ninact + nact + nsec
           do j = i, ninact + nact + nsec
               if (rank == 0) then
                   f(i, j) = DCMPLX(oner(i, j), onei(i, j))
               end if
               do k = 1, ninact + nelec

                   f(i, j) = f(i, j) + DCMPLX(int2r_f1(i, j, k, k), int2i_f1(i, j, k, k))
                   f(i, j) = f(i, j) - DCMPLX(int2r_f2(i, k, k, j), int2i_f2(i, k, k, j))

!Iwamuro modify
!                     write(*,*)f(i,j)
               End do           ! k

               f(j, i) = DCONJG(f(i, j))

           End do       ! j
       End do          ! i

       if (rank == 0) then ! Process limits for output
           write (3000, *) ' '
           write (3000, *) 'OFF DIAGONAL ELEMENTS OF FOCK MATRIX WHICH IS LARGER THAN 1.0d-06 '
           write (3000, *) ' '
           do i = 1, ninact + nact + nsec
           do j = i, ninact + nact + nsec
               if ((i /= j) .and. (ABS(f(i, j)) > 1.0d-6)) then
!            if(i/=j)then
                   write (3000, '(2I4,2E20.10)') i, j, f(i, j)
               end if
           end do
           end do
           write (3000, *) ' '
           write (3000, *) 'THESE DIAGONAL ELEMENTS SHOULD BE CORESPOND TO HF SPINOR ENERGY '
           write (3000, *) ' '
           write (3000, *) '  NO.   Spinor Energy(Re)   Spinor Energy(Im) '&
           &, 'Spinor Energy (HF)        ERROR'
           do i = 1, ninact + nact + nsec
               write (3000, '(I4,4E20.10)') i, f(i, i), orbmo(i), orbmo(i) - dble(f(i, i))
           end do

           write (3000, *) 'fockhf end'
       end if
   end SUBROUTINE fockhf1_ty
