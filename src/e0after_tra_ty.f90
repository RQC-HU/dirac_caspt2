! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE e0aftertra_ty

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_file_manager
    use four_caspt2_module

    Implicit NONE

    integer :: ii, jj, kk, ll
    integer :: j, i, k, l
    integer :: e0after_unit

    real*8 :: dr, di
    complex*16 :: oneeff, cmplxint, dens, energyHF(2)
    complex*16, allocatable :: energy(:, :)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    print *, 'EIGEN(1)', eigen(1)

    Allocate (energy(nroot, 4))
    energy(1:nroot, 1:4) = 0.0d+00

    debug = .FALSE.
    thres = 1.0d-15
    e0after_unit = default_unit
!        thres = 0.0d+00
    if (rank == 0) then
        call open_unformatted_file(unit=e0after_unit, file='e0after', status='new', optional_action='write')
!        AT PRESENT, CODE OF COMPLEX TYPE EXISTS !

        print *, 'iroot = ', iroot
    end if

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

        Call tramo1_ty(i, i, cmplxint)
!            write(*,'(I4,E20.10)')i,DBLE(cmplxint)
        energyHF(1) = energyHF(1) + cmplxint

    end do

!         print *,'energyHF(1)',energyHF(1)

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

            Call tramo2_ty(i, i, j, j, cmplxint)

            energyHF(2) = energyHF(2) + (0.5d+00)*cmplxint

            Call tramo2_ty(i, j, j, i, cmplxint)

            energyHF(2) = energyHF(2) - (0.5d+00)*cmplxint

        end do
    end do

    energyHF(2) = energyHF(2) + CONJG(energyHF(2))

!         print *,'energyHF(2)',energyHF(2)

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

        Call tramo1_ty(i, i, cmplxint)

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

            Call tramo2_ty(i, i, j, j, cmplxint)

            energy(iroot, 2) = energy(iroot, 2) + (0.5d+00)*cmplxint

            Call tramo2_ty(i, j, j, i, cmplxint)

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

                Call tramo2_ty(i, j, k, k, cmplxint)

                oneeff = oneeff + cmplxint

                Call tramo2_ty(i, k, k, j, cmplxint)

                oneeff = oneeff - cmplxint

300         end do           ! k

            Call tramo1_ty(i, j, cmplxint)

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
!             debug = .TRUE. ; print *, i,j,k,l
!          else
!             debug = .FALSE.
!          endif

                    Call tramo2_ty(i, j, k, l, cmplxint)

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
                        ! Only master rank are allowed to create files used by CASPT2 except for MDCINTNEW.
                        if (iroot == 1 .and. rank == 0) write (e0after_unit) i, j, k, l, DBLE(cmplxint), DBLE(dens)

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

100             end do        ! l
            end do    ! k
        end do       ! j
    end do          ! i

    energy(iroot, 4) = energy(iroot, 4) + CONJG(energy(iroot, 4))

!         if(ABS(eigen(iroot)-ecore &
!         -(energy(iroot,1)+energy(iroot,2)+energy(iroot,3)+energy(iroot,4))) &
!          > 1.0d-5 ) then
    if (rank == 0) then
        print *, 'energy 1 =', energy(iroot, 1)
        print *, 'energy 2 =', energy(iroot, 2)
        print *, 'energy 3 =', energy(iroot, 3)
        print *, 'energy 4 =', energy(iroot, 4)

        print *, iroot, 't-energy(1-4)', &
            energy(iroot, 1) + energy(iroot, 2) + energy(iroot, 3) + energy(iroot, 4)

        print *, iroot, 't-energy', &
            eigen(iroot) - ecore
        print *, iroot, 'eigen e0', &
            eigen(iroot)

        print *, 'C the error ', &
            eigen(iroot) - ecore &
            - (energy(iroot, 1) + energy(iroot, 2) + energy(iroot, 3) + energy(iroot, 4))
!         else
!            print *,'C the error ', &
!            eigen(iroot)-ecore &
!            -(energy(iroot,1)+energy(iroot,2)+energy(iroot,3)+energy(iroot,4))
!         end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      end do                    ! iroot = 1, nroot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        print *, 'energy HF  =', energyHF(1) + energyHF(2) + ecore
    end if

!!###   end do ! about type
    if (rank == 0) then  ! Only master ranks are allowed to create files used by CASPT2 except for MDCINTNEW.
        close (e0after_unit)
    end if
1000 continue
    deallocate (energy)
    print *, 'e0aftertra end'
End subroutine e0aftertra_ty

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE e0aftertrac_ty

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_file_manager
    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include "mpif.h"
#endif
    integer :: ii, jj, kk, ll
    integer :: j, i, k, l
    integer :: e0after_unit

    real*8 :: dr, di
    complex*16 :: oneeff, cmplxint, dens, energyHF(2)
    complex*16, allocatable :: energy(:, :)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    Allocate (energy(nroot, 4))
    energy(1:nroot, 1:4) = 0.0d+00

    debug = .FALSE.
    thres = 1.0d-15
    e0after_unit = default_unit
!        thres = 0.0d+00
    if (rank == 0) then ! Only master ranks are allowed to create files used by CASPT2 except for MDCINTNEW.
        call open_unformatted_file(unit=e0after_unit,file='e0after',status='new',optional_action='write')
!        AT PRESENT, CODE OF COMPLEX TYPE EXISTS !
        print *, 'iroot = ', iroot
    end if

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

    do i = rank + 1, ninact + nelec, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        ! do i = 1, ninact + nelec

        cmplxint = 0.0d+00

        Call tramo1_ty(i, i, cmplxint)
!            write(*,'(I4,E20.10)')i,DBLE(cmplxint)
        energyHF(1) = energyHF(1) + cmplxint

    end do
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, energyHF(1), 1, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) print *, 'energyHF(1)', energyHF(1)
!         do i = 1, ninact
!
!            cmplxint = 0.0d+00
!
!            Call tramo1_ty ( i, i, cmplxint)
!            energyHF(1) = energyHF(1) + cmplxint
!
!         end do
!
!         print *,'energyHF(1)',energyHF(1)
!
!         do i = ninact+1, ninact+nelec
!
!            cmplxint = 0.0d+00
!
!            Call tramo1_ty ( i, i, cmplxint)
!            energyHF(1) = energyHF(1) + cmplxint
!
!         end do
!
!         print *,'energyHF(1)',energyHF(1)

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

    do i = rank + 1, ninact + nelec, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        ! do i = 1, ninact + nelec
        do j = i, ninact + nelec

            Call tramo2_ty(i, i, j, j, cmplxint)

            energyHF(2) = energyHF(2) + (0.5d+00)*cmplxint

            Call tramo2_ty(i, j, j, i, cmplxint)

            energyHF(2) = energyHF(2) - (0.5d+00)*cmplxint

        end do
    end do

    energyHF(2) = energyHF(2) + DCONJG(energyHF(2))
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, energyHF(2), 1, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) print *, 'energyHF(2)', energyHF(2)

!Iwamuro modify
    if (rank == 0) print *, 'Iwamuro modify'

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!         energy 1            !
!"""""""""""""""""""""""""""""!
!   One-electron sumation     !
!                             !
!   Inactive (core) part      !
!                             !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!"""""""""""""""""""""""""""""
    do i = rank + 1, ninact, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        ! do i = 1, ninact

        Call tramo1_ty(i, i, cmplxint)

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
    do i = rank + 1, ninact, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        ! do i = 1, ninact
        do j = i, ninact

            Call tramo2_ty(i, i, j, j, cmplxint)

            energy(iroot, 2) = energy(iroot, 2) + (0.5d+00)*cmplxint

            Call tramo2_ty(i, j, j, i, cmplxint)

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
    do i = rank + ninact + 1, ninact + nact, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        ! do i = ninact + 1, ninact + nact
        do j = i, ninact + nact

            oneeff = 0.0d+00

            do k = 1, ninact            ! kk is inactive spinor

                Call tramo2_ty(i, j, k, k, cmplxint)

                oneeff = oneeff + cmplxint

                Call tramo2_ty(i, k, k, j, cmplxint)

                oneeff = oneeff - cmplxint

300         end do           ! k

            Call tramo1_ty(i, j, cmplxint)

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

    do i = rank + ninact + 1, ninact + nact, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        ! do i = ninact + 1, ninact + nact
        do j = ninact + 1, ninact + nact
            do k = ninact + 1, ninact + nact
                do l = i, ninact + nact

!          if((i < ninact+3).and.(j < ninact+3).and.(k < ninact+3).and.(l < ninact+3)) then
!             debug = .TRUE. ; print *, i,j,k,l
!          else
!             debug = .FALSE.
!          endif

                    Call tramo2_ty(i, j, k, l, cmplxint)

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
                        if (iroot == 1 .and. rank == 0) write (e0after_unit) i, j, k, l, DBLE(cmplxint), DBLE(dens) ! Only master ranks are allowed to create files used by CASPT2 except for MDCINTNEW.

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

100             end do        ! l
            end do    ! k
        end do       ! j
    end do          ! i

    energy(iroot, 4) = energy(iroot, 4) + CONJG(energy(iroot, 4))

!         if(ABS(eigen(iroot)-ecore &
!         -(energy(iroot,1)+energy(iroot,2)+energy(iroot,3)+energy(iroot,4))) &
!          > 1.0d-5 ) then
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, energy(iroot, 1), 1, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, energy(iroot, 2), 1, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, energy(iroot, 3), 1, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, energy(iroot, 4), 1, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    if (rank == 0) then
        print *, 'energy 1 =', energy(iroot, 1)
        print *, 'energy 2 =', energy(iroot, 2)
        print *, 'energy 3 =', energy(iroot, 3)
        print *, 'energy 4 =', energy(iroot, 4)

        print *, iroot, 't-energy(1-4)', &
            energy(iroot, 1) + energy(iroot, 2) + energy(iroot, 3) + energy(iroot, 4)

        print *, iroot, 't-energy', &
            eigen(iroot) - ecore
        print *, iroot, 'eigen e0', &
            eigen(iroot)

        print *, 'C the error ', &
            eigen(iroot) - ecore &
            - (energy(iroot, 1) + energy(iroot, 2) + energy(iroot, 3) + energy(iroot, 4))

! Iwamuro modify
        print *, 'Iwamuro modify'

!         else
!            print *,'C the error ', &
!            eigen(iroot)-ecore &
!            -(energy(iroot,1)+energy(iroot,2)+energy(iroot,3)+energy(iroot,4))
!         end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      end do                    ! iroot = 1, nroot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        print *, 'CAUTION! HF energy may not be obtained correctly '
        print *, 'energy HF  =', energyHF(1) + energyHF(2) + ecore
    end if
!!###   end do ! about type
    if (rank == 0) then ! Only master ranks are allowed to create files used by CASPT2 except for MDCINTNEW.
        close (e0after_unit)
    end if
1000 continue
    deallocate (energy)
!      print *,'e0aftertrac end'
! Iwamuro modify
    if (rank == 0) print *, 'e0aftertrac_ty end'
End subroutine e0aftertrac_ty
