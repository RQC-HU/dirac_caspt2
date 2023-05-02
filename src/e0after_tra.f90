! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE e0aftertra

! Subroutine to check whether the HF energy remains unchanged after MO transformation.

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_file_manager
    use four_caspt2_module

    Implicit NONE

    integer :: ii, jj, kk, ll
    integer :: j, i, k, l
    integer :: unit_e0after

    real*8 :: dr, di
    complex*16 :: oneeff, cmplxint, dens, energyHF(2)
    complex*16, allocatable :: energy(:, :)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    print *, 'EIGEN(1)', eigen(1)

    Allocate (energy(nroot, 4))
    energy(1:nroot, 1:4) = 0.0d+00

    debug = .FALSE.
    if (rank == 0) then
        call open_unformatted_file(unit=unit_e0after, file='e0after', status='replace', optional_action='write')
!        AT PRESENT, CODE OF COMPLEX TYPE EXISTS !

        print *, 'iroot = ', iroot
    end if

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
        energyHF(1) = energyHF(1) + cmplxint
    end do

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

            end do

            Call tramo1(i, j, cmplxint)

            oneeff = oneeff + cmplxint

            if (i == j) oneeff = 0.5d+00*oneeff

            if (realcvec) then

                ii = i - ninact
                jj = j - ninact
                Call dim1_density_R(ii, jj, dr)

                energy(iroot, 3) = energy(iroot, 3) + oneeff*dr

            else
                ii = i - ninact
                jj = j - ninact
                Call dim1_density(ii, jj, dr, di)

                dens = DCMPLX(dr, di)
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

                        dens = DCMPLX(dr, di)

                        ! Only master rank are allowed to create files used by CASPT2 except for MDCINTNEW.
                        if (iroot == 1 .and. rank == 0) write (unit_e0after) i, j, k, l, DBLE(cmplxint), DBLE(dens)

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

                            dens = DCMPLX(dr, di)
                            energy(iroot, 4) = energy(iroot, 4) &
                                               - (0.5d+00)*dens*cmplxint
                        end if

                    end if

                end do
            end do
        end do
    end do

    energy(iroot, 4) = energy(iroot, 4) + CONJG(energy(iroot, 4))

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

        print *, 'energy HF  =', energyHF(1) + energyHF(2) + ecore
    end if

    if (rank == 0) then  ! Only master ranks are allowed to create files used by CASPT2 except for MDCINTNEW.
        close (unit_e0after)
    end if
    deallocate (energy)
    print *, 'e0aftertra end'
End subroutine e0aftertra

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE e0aftertrac

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_file_manager
    use four_caspt2_module
#ifdef HAVE_MPI
    use module_mpi
#endif
    Implicit NONE

    integer :: ii, jj, kk, ll
    integer :: j, i, k, l

    real*8 :: dr, di
    complex*16 :: oneeff, cmplxint, dens, energyHF(2)
    complex*16, allocatable :: energy(:, :)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    Allocate (energy(nroot, 4))
    energy(1:nroot, 1:4) = 0.0d+00

    debug = .FALSE.
    if (rank == 0) print *, 'iroot = ', iroot

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

        cmplxint = 0.0d+00

        Call tramo1(i, i, cmplxint)
        energyHF(1) = energyHF(1) + cmplxint

    end do
#ifdef HAVE_MPI
    call allreduce_wrapper(mat=energyHF(1))
#endif
    if (rank == 0) print *, 'energyHF(1)', energyHF(1)

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
        do j = i, ninact + nelec

            Call tramo2(i, i, j, j, cmplxint)

            energyHF(2) = energyHF(2) + (0.5d+00)*cmplxint

            Call tramo2(i, j, j, i, cmplxint)

            energyHF(2) = energyHF(2) - (0.5d+00)*cmplxint

        end do
    end do

    energyHF(2) = energyHF(2) + DCONJG(energyHF(2))
#ifdef HAVE_MPI
    call allreduce_wrapper(mat=energyHF(2))
#endif
    if (rank == 0) print *, 'energyHF(2)', energyHF(2)

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
    do i = rank + 1, ninact, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
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
    do i = rank + ninact + 1, ninact + nact, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        do j = i, ninact + nact

            oneeff = 0.0d+00

            do k = 1, ninact            ! kk is inactive spinor
                Call tramo2(i, j, k, k, cmplxint)
                oneeff = oneeff + cmplxint

                Call tramo2(i, k, k, j, cmplxint)
                oneeff = oneeff - cmplxint
            end do           ! k

            Call tramo1(i, j, cmplxint)
            oneeff = oneeff + cmplxint
            if (i == j) oneeff = 0.5d+00*oneeff

            if (realcvec) then

                ii = i - ninact
                jj = j - ninact
                Call dim1_density_R(ii, jj, dr)

                energy(iroot, 3) = energy(iroot, 3) + oneeff*dr

            else
                ii = i - ninact
                jj = j - ninact
                Call dim1_density(ii, jj, dr, di)

                dens = DCMPLX(dr, di)
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
        do j = ninact + 1, ninact + nact
            do k = ninact + 1, ninact + nact
                do l = i, ninact + nact

                    Call tramo2(i, j, k, l, cmplxint)

                    If (i == l) cmplxint = cmplxint*(0.5d+00)
                    ii = i - ninact
                    jj = j - ninact
                    kk = k - ninact
                    ll = l - ninact
                    if (realcvec) then

                        Call dim2_density_R(ii, jj, kk, ll, dr)

                        energy(iroot, 4) = energy(iroot, 4) + (0.5d+00)*dr*cmplxint
                    else

                        Call dim2_density(ii, jj, kk, ll, dr, di)

                        dens = DCMPLX(dr, di)
                        energy(iroot, 4) = energy(iroot, 4) + (0.5d+00)*dens*cmplxint
                    end if

                    if (j == k) then

                        dr = 0.0d+00
                        di = 0.0d+00

                        if (realcvec) then

                            ii = i - ninact
                            ll = l - ninact

                            Call dim1_density_R(ii, ll, dr)

                            energy(iroot, 4) = energy(iroot, 4) - (0.5d+00)*dr*cmplxint
                        else

                            ii = i - ninact
                            ll = l - ninact

                            Call dim1_density(ii, ll, dr, di)

                            dens = DCMPLX(dr, di)
                            energy(iroot, 4) = energy(iroot, 4) - (0.5d+00)*dens*cmplxint
                        end if

                    end if

                end do
            end do
        end do
    end do

    energy(iroot, 4) = energy(iroot, 4) + CONJG(energy(iroot, 4))

#ifdef HAVE_MPI
    call allreduce_wrapper(mat=energy(iroot, 1:4))
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

        print *, 'CAUTION! HF energy may not be obtained correctly '
        print *, 'energy HF  =', energyHF(1) + energyHF(2) + ecore
    end if
    deallocate (energy)
    if (rank == 0) print *, 'e0aftertrac end'
End subroutine e0aftertrac
