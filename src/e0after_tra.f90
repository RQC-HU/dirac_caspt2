! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE e0aftertra

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE

       integer :: ii, jj, kk, ll, typetype
       integer :: j0, j, i, k, l, i0, i1, nuniq
       integer :: k0, l0, nint
       logical :: test

       real*8 :: i2r, i2i, dr, di, nsign
       complex*16 :: oneeff, cmplxint, dens, energyHF(2)
       complex*16, allocatable :: energy(:, :)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       write (*, *) 'EIGEN(1)', eigen(1)

       Allocate (energy(nroot, 4))
       energy(1:nroot, 1:4) = 0.0d+00

       debug = .FALSE.
       thres = 1.0d-15
!        thres = 0.0d+00

       open (5, file='e0after', status='unknown', form='unformatted')

!        AT PRESENT, CODE OF COMPLEX TYPE EXISTS !

       write (*, *) 'iroot = ', iroot

!        Do iroot = 1, nroot

!        Do iroot = 1, 1

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!         energy HF1          !
!"""""""""""""""""""""""""""""!
!   One-electron sumation     !
!                             !
!   Inactive (core) part      !
!                             !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!"""""""""""""""""""""""""""""
       energyHF(1) = 0.0d+00

       do i = 1, ninact + nelec

           cmplxint = 0.0d+00

           Call tramo1(i, i, cmplxint)
!            write(*,'(I4,E20.10)')i,DBLE(cmplxint)
           energyHF(1) = energyHF(1) + cmplxint

       end do

!         write(*,*)'energyHF(1)',energyHF(1)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!         energy HF2          !
!"""""""""""""""""""""""""""""!   1/2*[(rr|tt)-(rt|tr)}
!   Two-electron sumation     !
!                             !    for inactive r and t
!   Inactive (core) part      !
!                             !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!"""""""""""""""""""""""""""""
       energyHF(2) = 0.0d+00

       do i = 1, ninact + nelec
           do j = i, ninact + nelec

               Call tramo2(i, i, j, j, cmplxint)

               energyHF(2) = energyHF(2) + (0.5d+00)*cmplxint

               Call tramo2(i, j, j, i, cmplxint)

               energyHF(2) = energyHF(2) - (0.5d+00)*cmplxint

           end do
       end do

       energyHF(2) = energyHF(2) + CONJG(energyHF(2))

!         write(*,*)'energyHF(2)',energyHF(2)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!         energy 1            !
!"""""""""""""""""""""""""""""!
!   One-electron sumation     !
!                             !
!   Inactive (core) part      !
!                             !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!"""""""""""""""""""""""""""""
       do i = 1, ninact

           Call tramo1(i, i, cmplxint)

           energy(iroot, 1) = energy(iroot, 1) + cmplxint

       end do

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!         energy 2            !
!"""""""""""""""""""""""""""""!   1/2*[(rr|tt)-(rt|tr)}
!   Two-electron sumation     !
!                             !    for inactive r and t
!   Inactive (core) part      !
!                             !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!"""""""""""""""""""""""""""""
       do i = 1, ninact
           do j = i, ninact

               Call tramo2(i, i, j, j, cmplxint)

               energy(iroot, 2) = energy(iroot, 2) + (0.5d+00)*cmplxint

               Call tramo2(i, j, j, i, cmplxint)

               energy(iroot, 2) = energy(iroot, 2) - (0.5d+00)*cmplxint

           end do
       end do

       energy(iroot, 2) = energy(iroot, 2) + CONJG(energy(iroot, 2))

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!         energy 3            !
!"""""""""""""""""""""""""""""!
!   One-electron sumation     !
!                             !  hij + siguma [ (kk|ij)-(kj|ik) ]
!   Active part               !          k
!                             !
!   With effective one-e-int  !  hij + siguma [ (ij|kk)-(ik|kj) ]
!                             !          k
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!"""""""""""""""""""""""""""""
       do i = ninact + 1, ninact + nact
           do j = i, ninact + nact

               oneeff = 0.0d+00

               do k = 1, ninact            ! kk is inactive spinor

                   Call tramo2(i, j, k, k, cmplxint)

                   oneeff = oneeff + cmplxint

                   Call tramo2(i, k, k, j, cmplxint)

                   oneeff = oneeff - cmplxint

300            end do           ! k

               Call tramo1(i, j, cmplxint)

               oneeff = oneeff + cmplxint

!___________________________________________________________!
               !
               if (i == j) oneeff = 0.5d+00*oneeff             !
!___________________________________________________________!

               if (realcvec) then

                   ii = i - ninact
                   jj = j - ninact
                   Call dim1_density_R(ii, jj, dr)

                   energy(iroot, 3) = energy(iroot, 3) + oneeff*dr

               else
                   ii = i - ninact
                   jj = j - ninact
                   Call dim1_density(ii, jj, dr, di)

                   dens = CMPLX(dr, di, 16)
!                  write(*,'(2I4,2E20.10)') i, j,DBLE(oneeff), DBLE(dens)
                   energy(iroot, 3) = energy(iroot, 3) + oneeff*dens

               end if
           end do
       end do

       energy(iroot, 3) = energy(iroot, 3) + CONJG(energy(iroot, 3))

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!         energy 4            !
!"""""""""""""""""""""""""""""!   1/2*[(ij|kl)<0|EijEkl|0>-delta(jk)(ij|jl)<0|Eil|0>}
!   Two-electron sumation     !
!                             !   i,j,k and l are active spinors
!   active part               !
!                             !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!"""""""""""""""""""""""""""""

       do i = ninact + 1, ninact + nact
           do j = ninact + 1, ninact + nact
               do k = ninact + 1, ninact + nact
                   do l = i, ninact + nact

!          if((i < ninact+3).and.(j < ninact+3).and.(k < ninact+3).and.(l < ninact+3)) then
!             debug = .TRUE. ; write(*,*) i,j,k,l
!          else
!             debug = .FALSE.
!          endif

                       Call tramo2(i, j, k, l, cmplxint)

                       If (i == l) cmplxint = cmplxint*(0.5d+00)

                       if (realcvec) then
                           ii = i - ninact
                           jj = j - ninact
                           kk = k - ninact
                           ll = l - ninact

                           Call dim2_density_R(ii, jj, kk, ll, dr)

                           energy(iroot, 4) = energy(iroot, 4) &
                                              + (0.5d+00)*dr*cmplxint
                       else
                           ii = i - ninact
                           jj = j - ninact
                           kk = k - ninact
                           ll = l - ninact

                           Call dim2_density(ii, jj, kk, ll, dr, di)

                           dens = CMPLX(dr, di, 16)

!                  if(iroot==1) write(*,'(4I3,2E20.10)') i, j,k,l,DBLE(cmplxint), DBLE(dens)
                           if (iroot == 1) write (5) i, j, k, l, DBLE(cmplxint), DBLE(dens)

                           energy(iroot, 4) = energy(iroot, 4) &
                                              + (0.5d+00)*dens*cmplxint
                       end if

                       if (j == k) then

                           dr = 0.0d+00
                           di = 0.0d+00

                           if (realcvec) then

                               ii = i - ninact
                               ll = l - ninact

                               Call dim1_density_R(ii, ll, dr)

                               energy(iroot, 4) = energy(iroot, 4) &
                                                  - (0.5d+00)*dr*cmplxint
                           else

                               ii = i - ninact
                               ll = l - ninact

                               Call dim1_density(ii, ll, dr, di)

                               dens = CMPLX(dr, di, 16)
                               energy(iroot, 4) = energy(iroot, 4) &
                                                  - (0.5d+00)*dens*cmplxint
                           end if

                       end if

100                end do        ! l
               end do    ! k
           end do       ! j
       end do          ! i

       energy(iroot, 4) = energy(iroot, 4) + CONJG(energy(iroot, 4))

!         if(ABS(eigen(iroot)-ecore &
!         -(energy(iroot,1)+energy(iroot,2)+energy(iroot,3)+energy(iroot,4))) &
!          > 1.0d-5 ) then

       write (*, *) 'energy 1 =', energy(iroot, 1)
       write (*, *) 'energy 2 =', energy(iroot, 2)
       write (*, *) 'energy 3 =', energy(iroot, 3)
       write (*, *) 'energy 4 =', energy(iroot, 4)

       write (*, *) iroot, 't-energy(1-4)', &
           energy(iroot, 1) + energy(iroot, 2) + energy(iroot, 3) + energy(iroot, 4)

       write (*, *) iroot, 't-energy', &
           eigen(iroot) - ecore
       write (*, *) iroot, 'eigen e0', &
           eigen(iroot)

       write (*, *) 'C the error ', &
           eigen(iroot) - ecore &
           - (energy(iroot, 1) + energy(iroot, 2) + energy(iroot, 3) + energy(iroot, 4))

!         else
!            write(*,*)'C the error ', &
!            eigen(iroot)-ecore &
!            -(energy(iroot,1)+energy(iroot,2)+energy(iroot,3)+energy(iroot,4))
!         end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      end do                    ! iroot = 1, nroot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       write (*, *) 'energy HF  =', energyHF(1) + energyHF(2) + ecore

!!###   end do ! about type

       close (5)

1000   continue
       deallocate (energy)
       write (*, *) 'e0aftertra end'
   End subroutine e0aftertra

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE e0aftertrac

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE

       integer :: ii, jj, kk, ll, typetype
       integer :: j0, j, i, k, l, i0, i1, nuniq
       integer :: k0, l0, nint
       logical :: test

       real*8 :: i2r, i2i, dr, di, nsign
       complex*16 :: oneeff, cmplxint, dens, energyHF(2)
       complex*16, allocatable :: energy(:, :)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       Allocate (energy(nroot, 4))
       energy(1:nroot, 1:4) = 0.0d+00

       debug = .FALSE.
       thres = 1.0d-15
!        thres = 0.0d+00

       open (5, file='e0after', status='unknown', form='unformatted')

!        AT PRESENT, CODE OF COMPLEX TYPE EXISTS !

       write (*, *) 'iroot = ', iroot

!        Do iroot = 1, nroot

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!         energy HF1          !
!"""""""""""""""""""""""""""""!
!   One-electron sumation     !
!                             !
!   Inactive (core) part      !
!                             !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!"""""""""""""""""""""""""""""
       energyHF(1) = 0.0d+00

       do i = 1, ninact + nelec

           cmplxint = 0.0d+00

           Call tramo1(i, i, cmplxint)
!            write(*,'(I4,E20.10)')i,DBLE(cmplxint)
           energyHF(1) = energyHF(1) + cmplxint

       end do

!         do i = 1, ninact
!
!            cmplxint = 0.0d+00
!
!            Call tramo1 ( i, i, cmplxint)
!            energyHF(1) = energyHF(1) + cmplxint
!
!         end do
!
!         write(*,*)'energyHF(1)',energyHF(1)
!
!         do i = ninact+1, ninact+nelec
!
!            cmplxint = 0.0d+00
!
!            Call tramo1 ( i, i, cmplxint)
!            energyHF(1) = energyHF(1) + cmplxint
!
!         end do
!
!         write(*,*)'energyHF(1)',energyHF(1)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!         energy HF2          !
!"""""""""""""""""""""""""""""!   1/2*[(rr|tt)-(rt|tr)}
!   Two-electron sumation     !
!                             !    for inactive r and t
!   Inactive (core) part      !
!                             !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!"""""""""""""""""""""""""""""
       energyHF(2) = 0.0d+00

       do i = 1, ninact + nelec
           do j = i, ninact + nelec

               Call tramo2(i, i, j, j, cmplxint)
!               write(*,*)"tramo2 1"

               energyHF(2) = energyHF(2) + (0.5d+00)*cmplxint

               Call tramo2(i, j, j, i, cmplxint)
!               write(*,*)"tramo2 2"

               energyHF(2) = energyHF(2) - (0.5d+00)*cmplxint

           end do
       end do

       energyHF(2) = energyHF(2) + DCONJG(energyHF(2))

       write (*, *) 'energyHF(2)', energyHF(2)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!         energy 1            !
!"""""""""""""""""""""""""""""!
!   One-electron sumation     !
!                             !
!   Inactive (core) part      !
!                             !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!"""""""""""""""""""""""""""""
       do i = 1, ninact

           Call tramo1(i, i, cmplxint)

           energy(iroot, 1) = energy(iroot, 1) + cmplxint

       end do

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!         energy 2            !
!"""""""""""""""""""""""""""""!   1/2*[(rr|tt)-(rt|tr)}
!   Two-electron sumation     !
!                             !    for inactive r and t
!   Inactive (core) part      !
!                             !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!"""""""""""""""""""""""""""""
       do i = 1, ninact
           do j = i, ninact

               Call tramo2(i, i, j, j, cmplxint)
!               write(*,*)"tramo2 3"
               energy(iroot, 2) = energy(iroot, 2) + (0.5d+00)*cmplxint

               Call tramo2(i, j, j, i, cmplxint)
!               write(*,*)"tramo2 4"
               energy(iroot, 2) = energy(iroot, 2) - (0.5d+00)*cmplxint

           end do
       end do

       energy(iroot, 2) = energy(iroot, 2) + CONJG(energy(iroot, 2))

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!         energy 3            !
!"""""""""""""""""""""""""""""!
!   One-electron sumation     !
!                             !  hij + siguma [ (kk|ij)-(kj|ik) ]
!   Active part               !          k
!                             !
!   With effective one-e-int  !  hij + siguma [ (ij|kk)-(ik|kj) ]
!                             !          k
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!"""""""""""""""""""""""""""""
       do i = ninact + 1, ninact + nact
           do j = i, ninact + nact

               oneeff = 0.0d+00

               do k = 1, ninact            ! kk is inactive spinor

                   Call tramo2(i, j, k, k, cmplxint)
!                  write(*,*) "tramo2 5"
                   oneeff = oneeff + cmplxint

                   Call tramo2(i, k, k, j, cmplxint)
!                  write(*,*)"tramo2 6"
                   oneeff = oneeff - cmplxint

300            end do           ! k

               Call tramo1(i, j, cmplxint)

               oneeff = oneeff + cmplxint

!___________________________________________________________!
               !
               if (i == j) oneeff = 0.5d+00*oneeff             !
!___________________________________________________________!

               if (realcvec) then

                   ii = i - ninact
                   jj = j - ninact
                   Call dim1_density_R(ii, jj, dr)

                   energy(iroot, 3) = energy(iroot, 3) + oneeff*dr

               else
                   ii = i - ninact
                   jj = j - ninact
                   Call dim1_density(ii, jj, dr, di)

                   dens = CMPLX(dr, di, 16)
!                  write(*,'(2I4,2E20.10)') i, j,DBLE(oneeff), DBLE(dens)
                   energy(iroot, 3) = energy(iroot, 3) + oneeff*dens

               end if
           end do
       end do

       energy(iroot, 3) = energy(iroot, 3) + CONJG(energy(iroot, 3))
!Iwamuro modify
!         write(*,*)"energy(iroot,3)",energy(iroot,3)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!         energy 4            !
!"""""""""""""""""""""""""""""!   1/2*[(ij|kl)<0|EijEkl|0>-delta(jk)(ij|jl)<0|Eil|0>}
!   Two-electron sumation     !
!                             !   i,j,k and l are active spinors
!   active part               !
!                             !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!"""""""""""""""""""""""""""""

       do i = ninact + 1, ninact + nact
           do j = ninact + 1, ninact + nact
               do k = ninact + 1, ninact + nact
                   do l = i, ninact + nact

!          if((i < ninact+3).and.(j < ninact+3).and.(k < ninact+3).and.(l < ninact+3)) then
!             debug = .TRUE. ; write(*,*) i,j,k,l
!          else
!             debug = .FALSE.
!          endif

                       Call tramo2(i, j, k, l, cmplxint)
!                     write(*,*)"tramo2 7"
!Iwamuro modify
!                       write(*,*)'i, j, k, l, cmplxint ='
!                       write(*,'("testint2",4I4,2E15.5)') i, j, k, l, cmplxint

                       If (i == l) cmplxint = cmplxint*(0.5d+00)

                       if (realcvec) then
                           ii = i - ninact
                           jj = j - ninact
                           kk = k - ninact
                           ll = l - ninact

                           Call dim2_density_R(ii, jj, kk, ll, dr)
! Iwamuro modify
!                             write(*,*)'i, jj, kk, ll, dr ='
!                             write(*, '(4I4, E15.5)') ii, jj, kk, ll, dr

                           energy(iroot, 4) = energy(iroot, 4) &
                                              + (0.5d+00)*dr*cmplxint
                       else
                           ii = i - ninact
                           jj = j - ninact
                           kk = k - ninact
                           ll = l - ninact

                           Call dim2_density(ii, jj, kk, ll, dr, di)
! Iwamuro modify
!                             write(*,*)'ii, jj, kk, ll, dr, di ='
!                             Write(*,'(4i4, 2e15.5)') ii, jj, kk, ll, dr, di

                           dens = CMPLX(dr, di, 16)

!                  if(iroot==1) write(*,'(4I3,2E20.10)') i, j,k,l,DBLE(cmplxint), DBLE(dens)
                           if (iroot == 1) write (5) i, j, k, l, DBLE(cmplxint), DBLE(dens)

                           energy(iroot, 4) = energy(iroot, 4) &
                                              + (0.5d+00)*dens*cmplxint

!Iwamuro modify
!                           write(*,*) "energy(iroot,4)1", energy(iroot,4)
                       end if

                       if (j == k) then

                           dr = 0.0d+00
                           di = 0.0d+00

                           if (realcvec) then

                               ii = i - ninact
                               ll = l - ninact

                               Call dim1_density_R(ii, ll, dr)
! Iwamuro modify
!                             write(*,*)'i, ll, dr ='
!                             write(*,'(2I4, E15.5)') ii, ll, dr

                               energy(iroot, 4) = energy(iroot, 4) &
                                                  - (0.5d+00)*dr*cmplxint
!Iwamuro modify
!                           write(*,*) "energy(iroot,4)2", energy(iroot,4)
                           else

                               ii = i - ninact
                               ll = l - ninact

                               Call dim1_density(ii, ll, dr, di)
!Iwamuro modify
!                             write(*,*)'ii, ll, dr, di ='
!                             write(*, '(2I4, 2E15.5)') ii, ll, dr, di

                               dens = CMPLX(dr, di, 16)
                               energy(iroot, 4) = energy(iroot, 4) &
                                                  - (0.5d+00)*dens*cmplxint
!Iwamuro modify
!                           write(*,*) "energy(iroot,4)3", energy(iroot,4)
                           end if

                       end if

100                end do        ! l
               end do    ! k
           end do       ! j
       end do          ! i

       energy(iroot, 4) = energy(iroot, 4) + CONJG(energy(iroot, 4))

!         if(ABS(eigen(iroot)-ecore &
!         -(energy(iroot,1)+energy(iroot,2)+energy(iroot,3)+energy(iroot,4))) &
!          > 1.0d-5 ) then

       write (*, *) 'energy 1 =', energy(iroot, 1)
       write (*, *) 'energy 2 =', energy(iroot, 2)
       write (*, *) 'energy 3 =', energy(iroot, 3)
       write (*, *) 'energy 4 =', energy(iroot, 4)

       write (*, *) iroot, 't-energy(1-4)', &
           energy(iroot, 1) + energy(iroot, 2) + energy(iroot, 3) + energy(iroot, 4)

       write (*, *) iroot, 't-energy', &
           eigen(iroot) - ecore
       write (*, *) iroot, 'eigen e0', &
           eigen(iroot)

       write (*, *) 'C the error ', &
           eigen(iroot) - ecore &
           - (energy(iroot, 1) + energy(iroot, 2) + energy(iroot, 3) + energy(iroot, 4))

!         else
!            write(*,*)'C the error ', &
!            eigen(iroot)-ecore &
!            -(energy(iroot,1)+energy(iroot,2)+energy(iroot,3)+energy(iroot,4))
!         end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      end do                    ! iroot = 1, nroot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       write (*, *) 'CAUTION! HF energy may not be obtained correctly '
       write (*, *) 'energy HF  =', energyHF(1) + energyHF(2) + ecore

!!###   end do ! about type

       close (5)

1000   continue
       deallocate (energy)
       write (*, *) 'e0aftertrac end'
   End subroutine e0aftertrac