! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE fockcasci ! TO MAKE FOCK MATRIX for CASCI state

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

!! NOW MAKE FOCK MATRIX FOR CASCI STATE
!! fij = hij + SIGUMA_kl[<0|Ekl|0>{(ij|kl)-(il|kj)}

       f = 0.0d+00

       write (*, *) 'enter building fock matrix'

       do i = 1, ninact + nact
           do j = i, ninact + nact

               f(i, j) = DCMPLX(oner(i, j), onei(i, j))

               do k = 1, ninact

                   Call intmo2(i, j, k, k, cmplxint)
!Iwamuro modify
!                  If (i==4 .and. j==6.and. abs(cmplxint) > 1.0d-05) then
!                   write(*,'(4I4, 4E20.10)') i,j,k,k,cmplxint,f(i,j)
!                  endif
                   f(i, j) = f(i, j) + cmplxint

                   Call intmo2(i, k, k, j, cmplxint)
!Iwamuro modify
!                  If (i==4 .and. j==6.and. abs(cmplxint) > 1.0d-05) then
!                    write(*,'(4I4,4E20.10)') i,k,k,j,cmplxint,f(i,j)
!                  endif
                   f(i, j) = f(i, j) - cmplxint

               End do           ! k

               do k = ninact + 1, ninact + nact              ! ACTIVE SPACE
                   do l = ninact + 1, ninact + nact           ! ACTIVE SPACE

                       If (realcvec) then
                           Call dim1_density_R(k - ninact, l - ninact, dr)
                           Call intmo2(i, j, k, l, cmplxint)
                           f(i, j) = f(i, j) + dr*cmplxint
                           Call intmo2(i, l, k, j, cmplxint)
                           f(i, j) = f(i, j) - dr*cmplxint

                       Else
                           dr = 0.0d+00
                           Call dim1_density(k - ninact, l - ninact, dr, di)
                           dens = CMPLX(dr, di, 16)
!Iwamuro modify
!                        If (i==4 .and. j==6.and. abs(dens) > 0.0d-05)  Write( *, '("dens1",4I4, 2E20.10)') i,j,k,l,dens
!                        If (i==6 .and. j==4.and. abs(dens) > 0.0d-05)  Write( *, '("dens2",4I4, 2E20.10)') i,j,k,l,dens
                           Call intmo2(i, j, k, l, cmplxint)
                           f(i, j) = f(i, j) + dens*cmplxint
                           Call intmo2(i, l, k, j, cmplxint)
                           f(i, j) = f(i, j) - dens*cmplxint

                       End if

                   End do     ! l
               End do        ! k

               f(j, i) = DCONJG(f(i, j))
           end do       ! j
       end do          ! i

       do i = ninact + nact + 1, ninact + nact + nsec
           do j = i, ninact + nact + nsec

               f(i, j) = DCMPLX(oner(i, j), onei(i, j))
!               if(i==19.and.j==19)write(*,'("int1 ",2I4,2E20.10)')i,j,f(i,j)

               do k = 1, ninact

                   f(i, j) = f(i, j) + DCMPLX(int2r_f1(i, j, k, k), int2i_f1(i, j, k, k))
                   f(i, j) = f(i, j) - DCMPLX(int2r_f2(i, k, k, j), int2i_f2(i, k, k, j))

!                  if(i==19.and.j==19) write(*,'("+int2 ",4I4,2E20.10)') i,j,k,k, &
!                  & DCMPLX(int2r_f1(i,j,k,k),int2i_f1(i,j,k,k))
!
!                  if(i==19.and.j==19) write(*,'("-int2 ",4I4,2E20.10)') i,k,k,j, &
!                  & DCMPLX(int2r_f1(i,k,k,j),int2i_f1(i,k,k,j))

               End do           ! k

               do k = ninact + 1, ninact + nact              ! ACTIVE SPACE
                   do l = ninact + 1, ninact + nact           ! ACTIVE SPACE

                       If (realcvec) then
                           Call dim1_density_R(k - ninact, l - ninact, dr)

                           f(i, j) = f(i, j) + dr*DCMPLX(int2r_f1(i, j, k, l), int2i_f1(i, j, k, l))
                           f(i, j) = f(i, j) - dr*DCMPLX(int2r_f2(i, l, k, j), int2i_f2(i, l, k, j))

                       Else
                           Call dim1_density(k - ninact, l - ninact, dr, di)
                           dens = CMPLX(dr, di, 16)
!Iwamuro modify
!                        If (i==4 .and. j==6.and. abs(dens) > 0.0d-05)  Write( *, '("dens3",4I4, 2E20.10)') i,j,k,l,dens
!                        If (i==6 .and. j==4.and. abs(dens) > 0.0d-05)  Write( *, '("dens4",4I4, 2E20.10)') i,j,k,l,dens
                           f(i, j) = f(i, j) + dens*DCMPLX(int2r_f1(i, j, k, l), int2i_f1(i, j, k, l))
                           f(i, j) = f(i, j) - dens*DCMPLX(int2r_f2(i, l, k, j), int2i_f2(i, l, k, j))

!                        if(i==19.and.j==19) write(*,'("+int2 ",4I4,2E20.10)') i,j,k,l, &
!                        & DCMPLX(int2r_f1(i,j,k,l),int2i_f1(i,j,k,l))

!                        if(i==19.and.j==19) write(*,'("-int2 ",4I4,2E20.10)') i,l,k,j, &
!                        & DCMPLX(int2r_f2(i,l,k,j),int2i_f2(i,l,k,j))

!                        if(i==19.and.j==19) write(*,'("dens ",2I4,2E20.10)') k,l, dens

                       End if

                   End do     ! l
               End do        ! k

!Iwamuro modify
               !                 If (i==4 .and. j==6)  Write( *, '(4I4, 2E20.10)') i,j,k,l,f(i,j)
               f(j, i) = DCONJG(f(i, j))
!Iwamuro modify
!                  write(*,'("fock",2I4,2E20.10)')i,j,f(j,i)
!                  if(i==19.and.j==19)write(*,'("fock ",2I4,2E20.10)')i,j,f(i,j)
!                  write(*,'("fock",2I4,2E20.10)')i,j,f(i,j)

           end do       ! j
       end do          ! i

       write (*, *) 'fockcasci end'
   end
