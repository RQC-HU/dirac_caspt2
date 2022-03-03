! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE fockcasci_ty ! TO MAKE FOCK MATRIX for CASCI state

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE
#ifdef HAVE_MPI
       include 'mpif.h'
#endif
       integer :: ii, jj, kk, ll
       integer :: j, i, k, l
       integer :: nint, n

       real*8 :: i2r, i2i, dr, di, nsign
       complex*16 :: cmplxint, dens
       integer :: datetmp0, datetmp1
       real(8) :: tsectmp0, tsectmp1
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!! NOW MAKE FOCK MATRIX FOR CASCI STATE
!! fij = hij + SIGUMA_kl[<0|Ekl|0>{(ij|kl)-(il|kj)}
       if (rank == 0) then ! Process limits for output
           write (*, *) 'enter building fock matrix'
       end if
       datetmp0 = initdate
       tsectmp0 = inittime
       if (rank == 0) call timing(datetmp0, tsectmp0, datetmp1, tsectmp1)
       datetmp0 = datetmp1
       tsectmp0 = tsectmp1

       f = 0.0d+00

       if (rank == 0) then ! Process limits for output
           write (*, *) 'enter building fock matrix'
       end if
       !$OMP parallel private(i,j,k,l,cmplxint,dr,di,dens)
       !$OMP do schedule(dynamic,2)
       do i = rank + 1, ninact + nact, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
           do j = i, ninact + nact
               f(i, j) = DCMPLX(oner(i, j), onei(i, j))
               do k = 1, ninact

                   !    Call intmo2_ty(i, j, k, k, cmplxint)

                   !    f(i, j) = f(i, j) + cmplxint
                   f(i, j) = f(i, j) + CMPLX(inttwr(i, j, k, k), inttwi(i, j, k, k), 16)
                   !    Call intmo2_ty(i, k, k, j, cmplxint)

                   !    f(i, j) = f(i, j) - cmplxint
                   f(i, j) = f(i, j) - CMPLX(inttwr(i, k, k, j), inttwi(i, k, k, j), 16)

               End do           ! k
               do k = ninact + 1, ninact + nact              ! ACTIVE SPACE
                   do l = ninact + 1, ninact + nact           ! ACTIVE SPACE

                       If (realcvec) then
                           Call dim1_density_R(k - ninact, l - ninact, dr)
                           !    Call intmo2_ty(i, j, k, l, cmplxint)
                           !    f(i, j) = f(i, j) + dr*cmplxint
                           f(i, j) = f(i, j) + dr*CMPLX(inttwr(i, j, k, l), inttwi(i, j, k, l), 16)
                           !    Call intmo2_ty(i, l, k, j, cmplxint)
                           !    f(i, j) = f(i, j) - dr*cmplxint
                           f(i, j) = f(i, j) - dr*CMPLX(inttwr(i, l, k, j), inttwi(i, l, k, j), 16)
                       Else
                           dr = 0.0d+00
                           Call dim1_density(k - ninact, l - ninact, dr, di)

                           dens = CMPLX(dr, di, 16)
                           !    Call intmo2_ty(i, j, k, l, cmplxint)
                           !    f(i, j) = f(i, j) + dens*cmplxint
                           f(i, j) = f(i, j) + dens*CMPLX(inttwr(i, j, k, l), inttwi(i, j, k, l), 16)

                           !    Call intmo2_ty(i, l, k, j, cmplxint)
                           !    f(i, j) = f(i, j) - dens*cmplxint
                           f(i, j) = f(i, j) - dens*CMPLX(inttwr(i, l, k, j), inttwi(i, l, k, j), 16)

                       End if

                   End do     ! l
               End do        ! k

               f(j, i) = DCONJG(f(i, j))
           end do       ! j
       end do          ! i
       !$OMP end do

       !$OMP do schedule(dynamic,2)
       do i = ninact + nact + 1 + rank, ninact + nact + nsec, nprocs
           do j = i, ninact + nact + nsec
               !    if (rank == 0) then
               f(i, j) = DCMPLX(oner(i, j), onei(i, j))
!               if(i==19.and.j==19)write(*,'("int1 ",2I4,2E20.10)')i,j,f(i,j)
               !    end if
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

               f(j, i) = DCONJG(f(i, j))
!                  if(i==19.and.j==19)write(*,'("fock ",2I4,2E20.10)')i,j,f(i,j)
!                  write(*,'("fock",2I4,2E20.10)')i,j,f(i,j)

           end do       ! j
       end do          ! i
       !$OMP end do
       !$OMP end parallel
#ifdef HAVE_MPI
       call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
       if (rank == 0) then  ! Process limits for output
           write (*, *) 'fockcasci before f allreduce'
       end if
       if (rank == 0) call timing(datetmp0, tsectmp0, datetmp1, tsectmp1)
       datetmp0 = datetmp1
       tsectmp0 = tsectmp1
#ifdef HAVE_MPI
       call MPI_Allreduce(MPI_IN_PLACE, f(1, 1), nmo**2, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
       if (rank == 0) then ! Process limits for output
           write (*, *) 'fockcasci end'
       end if
   end
