
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE e0aftertra

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_file_manager
    use module_global_variables
    use module_index_utils, only: convert_global_to_active_idx
    use module_realonly, only: realonly
#ifdef HAVE_MPI
    use module_mpi
#endif
    Implicit NONE

    integer :: ii, jj, kk, ll
    integer :: j, i, k, l

    real(8) :: dr, di
    complex*16 :: oneeff, cmplxint, dens, energyHF(2)
    complex*16, allocatable :: energy(:, :)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    Allocate (energy(nroot, 4))
    ! Initialize variables
    energyHF(:) = 0.0d+00
    energy(:, :) = 0.0d+00
    dr = 0.0d+00
    di = 0.0d+00
    oneeff = 0.0d+00
    cmplxint = 0.0d+00
    dens = 0.0d+00
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
    do i = rank + 1, ninact + nelec, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
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
    do i = rank + 1, ninact + nelec, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        do j = i, ninact + nelec
            Call tramo2(i, i, j, j, cmplxint)
            energyHF(2) = energyHF(2) + cmplxint

            Call tramo2(i, j, j, i, cmplxint)
            energyHF(2) = energyHF(2) - cmplxint
        end do
    end do
    energyHF(2) = 0.5d+00*energyHF(2)
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
            energy(iroot, 2) = energy(iroot, 2) + cmplxint

            Call tramo2(i, j, j, i, cmplxint)
            energy(iroot, 2) = energy(iroot, 2) - cmplxint
        end do
    end do
    energy(iroot, 2) = 0.5d+00*energy(iroot, 2)
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

            ii = convert_global_to_active_idx(i)
            jj = convert_global_to_active_idx(j)
            if (realonly%is_realonly()) then
                Call dim1_density_R(ii, jj, dr) ! di is always zero
            else
                Call dim1_density(ii, jj, dr, di)
            end if
            dens = DCMPLX(dr, di)
            energy(iroot, 3) = energy(iroot, 3) + oneeff*dens
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

                    ii = convert_global_to_active_idx(i)
                    jj = convert_global_to_active_idx(j)
                    kk = convert_global_to_active_idx(k)
                    ll = convert_global_to_active_idx(l)
                    if (realonly%is_realonly()) then
                        Call dim2_density_R(ii, jj, kk, ll, dr) ! di is always zero
                    else
                        Call dim2_density(ii, jj, kk, ll, dr, di)
                    end if
                    dens = DCMPLX(dr, di)
                    energy(iroot, 4) = energy(iroot, 4) + dens*cmplxint

                    if (j == k) then
                        if (realonly%is_realonly()) then
                            Call dim1_density_R(ii, ll, dr) ! di is always zero
                        else
                            Call dim1_density(ii, ll, dr, di)
                        end if
                        dens = DCMPLX(dr, di)
                        energy(iroot, 4) = energy(iroot, 4) - dens*cmplxint
                    end if
                end do
            end do
        end do
    end do
    energy(iroot, 4) = 0.5d+00*energy(iroot, 4)
    energy(iroot, 4) = energy(iroot, 4) + CONJG(energy(iroot, 4))

#ifdef HAVE_MPI
    call allreduce_wrapper(mat=energyHF(:))
    call allreduce_wrapper(mat=energy(iroot, :))
#endif

    if (rank == 0) then
        print '(a,x,i0)', 'selectroot =', iroot
        print *, 'core energy =', ecore
        print *, 'energyHF(1)', energyHF(1)
        print *, 'energyHF(2)', energyHF(2)
        print *, 'energy 1 =', energy(iroot, 1)
        print *, 'energy 2 =', energy(iroot, 2)
        print *, 'energy 3 =', energy(iroot, 3)
        print *, 'energy 4 =', energy(iroot, 4)
        print *, iroot, 't-energy(1-4)', sum(energy(iroot, :))
        print *, iroot, 't-energy', eigen(iroot) - ecore
        print *, iroot, 'eigen e0', eigen(iroot)
        print *, 'C the error ', eigen(iroot) - ecore - sum(energy(iroot, :))
        print *, 'CAUTION! HF energy may not be obtained correctly '
        print *, 'energy HF  =', sum(energyHF) + ecore
    end if
    deallocate (energy)
    if (rank == 0) print *, 'e0aftertra end'
End subroutine e0aftertra