! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE e0test_v2 ! test to calculate <i|H|i>=Ei i is solution of the CASCI

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer :: ii, jj, kk, ll
    integer :: j, i, k, l
    integer :: nint

    real*8 :: i2r, i2i, dr, di, nsign
    complex*16 :: oneeff, cmplxint, dens, energyHF(2)
    complex*16, allocatable :: energy(:, :)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    if (rank == 0) then ! Process limits for output
        write (*, *) "enter e0test"
    end if
    Allocate (energy(nroot, 4)); Call memplus(KIND(energy), SIZE(energy), 1)
    energy(:, :) = 0.0d+00
    debug = .TRUE.

    if (realc) then
!~~~~~~~~~~~~~~~~~~~~~~

!RRRRRRRRRRRRRRRRRRRRRRRRRRRRR!
!         energy 1            !
!"""""""""""""""""""""""""""""!   1/2*[(rr|tt)-(rt|tr)}
!   One-electron sumation     !
!                             !
!   Inactive (core) part      !
!                             !
!RRRRRRRRRRRRRRRRRRRRRRRRRRRRR!
!"""""""""""""""""""""""""""""!

        do i = 1, ninact
            ii = i
            !! Adding one-electron integral to the fock matrics is executed only by the master process
            !! because DIRAC's one-electron integral file (MRCONEE) is not
            !! devided even if DIRAC is executed in parallel (MPI).
            if (rank == 0) then
                energy(iroot, 1) = energy(iroot, 1) + oner(ii, ii)
            end if
        end do

        if (rank == 0) then ! Process limits for output
            write (*, *) 'energy 1 =', energy(iroot, 1)
        end if

!RRRRRRRRRRRRRRRRRRRRRRRRRRRRR!
!         energy 2            !
!"""""""""""""""""""""""""""""!   1/2*[(rr|tt)-(rt|tr)}
!   Two-electron sumation     !
!                             !    for inactive r and t
!   Inactive (core) part      !
!                             !
!RRRRRRRRRRRRRRRRRRRRRRRRRRRRR!
!"""""""""""""""""""""""""""""
        do i = 1, ninact
            do j = 1, ninact
                ii = i
                jj = j

                i2r = inttwr(ii, ii, jj, jj)
                energy(iroot, 2) = energy(iroot, 2) + (0.5d+00)*i2r

                i2r = inttwr(ii, jj, jj, ii)
                energy(iroot, 2) = energy(iroot, 2) - (0.5d+00)*i2r

            end do
        end do

        if (rank == 0) then ! Process limits for output
            write (*, *) 'energy 2 =', energy(iroot, 2)
        end if

!RRRRRRRRRRRRRRRRRRRRRRRRRRRRR!
!         energy 3            !
!"""""""""""""""""""""""""""""!  hij + siguma [ (kk|ij)-(kj|ik) ]
!   One-electron sumation     !          k
!                             !
!   Active part               !
!                             !
!   With effective one-e-int  !
!                             !
!RRRRRRRRRRRRRRRRRRRRRRRRRRRRR!
!"""""""""""""""""""""""""""""

        do i = ninact + 1, ninact + nact
            do j = ninact + 1, ninact + nact
                ii = i
                jj = j

                oneeff = 0.0d+00

                do k = 1, ninact            ! kk is inactive spinor
                    kk = k

                    i2r = inttwr(kk, kk, ii, jj)
                    oneeff = oneeff + i2r

                    i2r = inttwr(kk, jj, ii, kk)
                    oneeff = oneeff - i2r

                end do ! kk

                if (realcvec) then
                    Call dim1_density_R(ii, jj, dr)
                    if (rank == 0) then
                        energy(iroot, 3) = energy(iroot, 3) + (oner(ii, jj) + oneeff)*dr
                    end if
                else
                    Call dim1_density(ii, jj, dr, di)
                    if (rank == 0) then
                        energy(iroot, 3) = energy(iroot, 3) + (oner(ii, jj) + oneeff)*DCMPLX(dr, di)
                    end if
                end if
            end do
        end do

        if (rank == 0) then ! Process limits for output
            write (*, *) 'energy 3 =', energy(iroot, 3)
        end if

!RRRRRRRRRRRRRRRRRRRRRRRRRRRRR!
!         energy 4            !
!"""""""""""""""""""""""""""""!   1/2*[(ij|kl)<0|EijEkl|0>-delta(jk)(ij|jl)<0|Eil|0>}
!   Two-electron sumation     !
!                             !   i,j,k and l are ative spinors
!   active part               !
!                             !
!RRRRRRRRRRRRRRRRRRRRRRRRRRRRR!
!"""""""""""""""""""""""""""""

        do i = ninact + 1, ninact + nact
            do j = ninact + 1, ninact + nact
                do k = ninact + 1, ninact + nact
                    do l = ninact + 1, ninact + nact
                        ii = i
                        jj = j
                        kk = k
                        ll = l

                        i2r = 0.0d+00
                        i2i = 0.0d+00
                        dr = 0.0d+00
                        di = 0.0d+00

                        i2r = inttwr(ii, jj, kk, ll)

                        if (realcvec) then
                            Call dim2_density_R(ii, jj, kk, ll, dr)
                            energy(iroot, 4) = energy(iroot, 4) + (0.5d+00)*dr*i2r
                        else
                            Call dim2_density(ii, jj, kk, ll, dr, di)
                            energy(iroot, 4) = energy(iroot, 4) + (0.5d+00)*DCMPLX(dr, di)*i2r
                        end if

                        if (jj == kk) then

                            dr = 0.0d+00
                            di = 0.0d+00

                            if (realcvec) then
                                Call dim1_density_R(ii, ll, dr)
                                energy(iroot, 4) = energy(iroot, 4) - (0.5d+00)*dr*i2r
                            else
                                Call dim1_density(ii, ll, dr, di)
                                energy(iroot, 4) = energy(iroot, 4) - (0.5d+00)*DCMPLX(dr, di)*i2r
                            end if

                        end if

                    end do ! ll
                end do    ! kk
            end do       ! jj
        end do          ! ii

        if (rank == 0) then ! Process limits for output
            write (*, *) 'energy 4 =', energy(iroot, 4)

            write (*, *) iroot, 't-energy(1-4)', &
                energy(iroot, 1) + energy(iroot, 2) + energy(iroot, 3) + energy(iroot, 4)

            write (*, *) iroot, 't-energy ', &
                eigen(iroot) - ecore

            write (*, *) 'R the error ', eigen(iroot) - ecore &
                - (energy(iroot, 1) + energy(iroot, 2) + energy(iroot, 3) + energy(iroot, 4))
        end if
    else

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
        ii = 0
        do i = rank + 1, ninact + nelec, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
            cmplxint = 0.0d+00
            cmplxint = CMPLX(oner(i, i), onei(i, i), 16)
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

        do i = rank + 1, ninact + nelec, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
            do j = i, ninact + nelec

                cmplxint = 0.0d+00
                i2r = 0.0d+00
                i2i = 0.0d+00
                nsign = 0.0d+00

                i2r = inttwr(i, i, j, j)
                i2i = inttwi(i, i, j, j)

                cmplxint = CMPLX(i2r, i2i, 16)
                energyHF(2) = energyHF(2) + (0.5d+00)*cmplxint

                cmplxint = 0.0d+00
                i2r = 0.0d+00
                i2i = 0.0d+00
                nsign = 0.0d+00
                nint = 0

                i2r = inttwr(i, j, j, i)
                i2i = inttwi(i, j, j, i)

                cmplxint = CMPLX(i2r, i2i, 16)
                energyHF(2) = energyHF(2) - (0.5d+00)*cmplxint

            end do
        end do

        energyHF(2) = energyHF(2) + DCONJG(energyHF(2))

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
            cmplxint = CMPLX(oner(i, i), onei(i, i), 16)
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
                i2r = inttwr(i, i, j, j)
                i2i = inttwi(i, i, j, j)

                cmplxint = CMPLX(i2r, i2i, 16)
                energy(iroot, 2) = energy(iroot, 2) + (0.5d+00)*cmplxint

                i2r = inttwr(i, j, j, i)
                i2i = inttwi(i, j, j, i)

                cmplxint = CMPLX(i2r, i2i, 16)
                energy(iroot, 2) = energy(iroot, 2) - (0.5d+00)*cmplxint

            end do
        end do

        energy(iroot, 2) = energy(iroot, 2) + DCONJG(energy(iroot, 2))

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

                    i2r = inttwr(i, j, k, k)
                    i2i = inttwi(i, j, k, k)

                    cmplxint = CMPLX(i2r, i2i, 16)
                    oneeff = oneeff + cmplxint

                    i2r = inttwr(i, k, k, j)
                    i2i = inttwi(i, k, k, j)

                    cmplxint = CMPLX(i2r, i2i, 16)
                    oneeff = oneeff - cmplxint

                end do           ! kk
                cmplxint = CMPLX(oner(i, j), onei(i, j), 16)
                oneeff = oneeff + cmplxint

                if (i == j) oneeff = oneeff*(0.5d+00)

                if (realcvec) then
                    ii = i - ninact
                    jj = j - ninact
                    Call dim1_density_R(ii, jj, dr)
                    ii = i
                    jj = j

                    energy(iroot, 3) = energy(iroot, 3) + oneeff*dr

                else

                    Call dim1_density(i - ninact, j - ninact, dr, di)

                    dens = CMPLX(dr, di, 16)
                    energy(iroot, 3) = energy(iroot, 3) + oneeff*dens

                end if
            end do
        end do

        energy(iroot, 3) = energy(iroot, 3) + DCONJG(energy(iroot, 3))

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
                    do l = ninact + 1, ninact + nact

                        i2r = 0.0d+00
                        i2i = 0.0d+00
                        dr = 0.0d+00
                        di = 0.0d+00

                        i2r = inttwr(i, j, k, l)
                        i2i = inttwi(i, j, k, l)

                        cmplxint = CMPLX(i2r, i2i, 16)

                        if (realcvec) then
                            ii = i - ninact
                            jj = j - ninact
                            kk = k - ninact
                            ll = l - ninact

                            Call dim2_density_R(ii, jj, kk, ll, dr)
                            ii = i
                            jj = j
                            kk = k
                            ll = l

                            energy(iroot, 4) = energy(iroot, 4) + (0.5d+00)*dr*cmplxint
                        else
                            ii = i - ninact
                            jj = j - ninact
                            kk = k - ninact
                            ll = l - ninact

                            Call dim2_density(ii, jj, kk, ll, dr, di)

                            dens = CMPLX(dr, di, 16)

                            energy(iroot, 4) = energy(iroot, 4) + (0.5d+00)*dens*cmplxint
                        end if

                        if (j == k) then

                            dr = 0.0d+00
                            di = 0.0d+00

                            if (realcvec) then

                                ii = i - ninact
                                ll = l - ninact

                                Call dim1_density_R(ii, ll, dr)
                                ii = i
                                jj = j
                                kk = k
                                ll = l

                                energy(iroot, 4) = energy(iroot, 4) - (0.5d+00)*dr*cmplxint
                            else

                                ii = i - ninact
                                ll = l - ninact

                                Call dim1_density(ii, ll, dr, di)

                                ii = i
                                jj = j
                                kk = k
                                ll = l

                                dens = CMPLX(dr, di, 16)
                                energy(iroot, 4) = energy(iroot, 4) - (0.5d+00)*dens*cmplxint
                            end if

                        end if

                    end do        ! ll
                end do    ! kk
            end do       ! jj
        end do          ! ii

#ifdef HAVE_MPI
        call MPI_Allreduce(MPI_IN_PLACE, energy(iroot, 1), 1, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, energy(iroot, 2), 1, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, energy(iroot, 3), 1, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, energy(iroot, 4), 1, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, energyHF(1), 1, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, energyHF(2), 1, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

        if (rank == 0) then ! Process limits for output
            write (*, *) 'energy 1 =', energy(iroot, 1)
            write (*, *) 'energy 2 =', energy(iroot, 2)
            write (*, *) 'energy 3 =', energy(iroot, 3)
            write (*, *) 'energy 4 =', energy(iroot, 4)

            write (*, *) iroot, 't-energy(1-4)', &
                energy(iroot, 1) + energy(iroot, 2) + energy(iroot, 3) + energy(iroot, 4)

            write (*, *) iroot, 't-energy', &
                eigen(iroot) - ecore

            write (*, *) 'C the error ', &
                eigen(iroot) - ecore &
                - (energy(iroot, 1) + energy(iroot, 2) + energy(iroot, 3) + energy(iroot, 4))

        end if
        if (rank == 0) then ! Process limits for output
            write (*, *) 'energy HF  =', energyHF(1) + energyHF(2) + ecore
        end if
    end if
    deallocate (energy); Call memminus(KIND(energy), SIZE(energy), 1)
    if (rank == 0) then ! Process limits for output
        write (*, *) 'e0test end'
    end if
End subroutine e0test_v2
