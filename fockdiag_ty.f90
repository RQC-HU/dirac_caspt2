! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE fockdiag_ty

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE

       integer                 ::  i, j
       integer                 :: i0, j0, n, dimn, n0, n1, nspace(3, 3)
       logical                 :: test, cutoff

       complex*16              :: trace1, trace2
       real*8, allocatable :: fa(:, :)
       complex*16, allocatable :: fac(:, :), readmo(:, :, :)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       if (rank == 0) then ! Process limits for output
           write (normaloutput, *) 'fockdiag start'
       end if
       REALF = .TRUE.

       Do i = 1, ninact + nact + nsec
           Do j = 1, ninact + nact + nsec
               If (ABS(DIMAG(f(i, j))) > 1.0d-12) then
                   REALF = .FALSE.
               End if
           End do
       End do

       REALF = .FALSE.
       if (rank == 0) then ! Process limits for output
           write (normaloutput, *) 'REALF', REALF
       end if

       If (REALF) then          ! real*8
           Allocate (fa(nmo, nmo)); Call memplus(KIND(fa), SIZE(fa), 1)
           eps = 0.0d+00
           fa = 0.0d+00
       Else
           Allocate (fac(nmo, nmo)); Call memplus(KIND(fac), SIZE(fac), 2)
           eps = 0.0d+00
           fac = 0.0d+00
       End if

       nspace(1, 1) = 1
       nspace(2, 1) = ninact
       nspace(3, 1) = ninact

       nspace(1, 2) = ninact + 1
       nspace(2, 2) = ninact + nact
       nspace(3, 2) = nact

       nspace(1, 3) = ninact + nact + 1
       nspace(2, 3) = ninact + nact + nsec
       nspace(3, 3) = nsec

       Do i0 = 1, 3

           n0 = nspace(1, i0)
           n1 = nspace(2, i0)
           n = nspace(3, i0)

           if (rank == 0) then ! Process limits for output
               if (i0 == 1) write (normaloutput, *) 'FOR INACTIVE-INACTIVE ROTATION !'
               if (i0 == 2) write (normaloutput, *) 'FOR ACTIVE-ACTIVE ROTATION !'
               if (i0 == 3) write (normaloutput, *) 'FOR SECONDARY-SECONDARY ROTATION !'
           end if
           if (REALF) then

               Call rdiag0(n, n0, n1, fa(n0:n1, n0:n1), eps(n0:n1))

               write (5) n0, n1, n
               write (5) fa(n0:n1, n0:n1)
               if (rank == 0) then ! Process limits for output
                   write (normaloutput, *) n0, n1, n

                   write (normaloutput, *) 'fa '

                   do i = n0, n1
                       write (normaloutput, '(30E13.5)') (fa(i, j), j=n0, n1)
                   end do

                   write (normaloutput, *) 'f '
                   do i = n0, n1
                       write (normaloutput, '(30E13.5)') (DBLE(f(i, j)), j=n0, n1)
                   end do
               end if
           else
               Call cdiag0(n, n0, n1, fac(n0:n1, n0:n1), eps(n0:n1))

           end if

       End do ! i0

       if (REALF) then

           Call traci(fa(ninact + 1:ninact + nact, ninact + 1:ninact + nact))

           f(1:nmo, 1:nmo) = fa(1:nmo, 1:nmo)

           Call e0aftertra_ty

           deallocate (fa); Call memminus(KIND(fa), SIZE(fa), 1)

       else

           Call tracic(fac(ninact + 1:ninact + nact, ninact + 1:ninact + nact))

           f(1:nmo, 1:nmo) = fac(1:nmo, 1:nmo)

           Call e0aftertrac_ty

           deallocate (fac); Call memminus(KIND(fac), SIZE(fac), 2)

       end if

!         Do i0 = (ninact+nact)/2+1, nmo/2
!         Do j0 = (ninact+nact)/2+1, nmo/2
!
!            if(ABS(f(2*i0,2*j0)-DCONJG(f(2*i0-1,2*j0-1))) > 1.0d-10) then
!               write(*,'(2I4,2E20.10)')2*i0,2*j0,f(2*i0,2*j0)
!               write(*,'(2I4,2E20.10)')2*i0-1,2*j0-1,f(2*i0-1,2*j0-1)
!               write(*,*)' '
!            Endif
!
!         Enddo
!         Enddo
       if (rank == 0) then ! Only master ranks are allowed to create files used by CASPT2 except for MDCINTNEW.
           open (5, file='TRANSFOCK', status='unknown', form='unformatted')
           write (5) nmo
           write (5) f(1:nmo, 1:nmo)
           close (5)
       end if

       goto 1000
10     if (rank == 0) then ! Process limits for output
           write (normaloutput, *) 'reading err in orbcoeff'
       end if
1000   continue
       if (rank == 0) then ! Process limits for output
           write (normaloutput, *) 'fockdiag end'
       end if
   end subroutine fockdiag_ty
