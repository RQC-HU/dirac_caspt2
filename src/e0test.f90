! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE e0test ! test to calculate <i|H|i>=Ei i is solution of the CASCI

! Recalculates the CASCI energy(ecore+energy(iroot,1:4))
! using 1,2 electron integrals and CI coefficients
! and compare the eigenvalue(CASCI energy, eigen(iroot))
! obtained by diagonalizing the CASCI matrix

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables
    use module_index_utils, only: convert_global_to_active_idx
    use module_realonly
#ifdef HAVE_MPI
    use module_mpi
#endif
    Implicit NONE

    integer :: ii, jj, kk, ll
    integer :: j, i, k, l
    real(8) :: i2r, i2i, dr, di
    complex*16 :: oneeff, cmplxint, dens, energyHF(2)
    complex*16, allocatable :: energy(:, :)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    if (rank == 0) print *, "enter e0test"
    Allocate (energy(nroot, 4)); Call memplus(KIND(energy), SIZE(energy), 1)

    ! Initialize variables
    energyHF(:) = 0.0d+00
    energy(:, :) = 0.0d+00
    debug = .TRUE.
    cmplxint = 0.0d+00
    i2r = 0.0d+00
    i2i = 0.0d+00
    dr = 0.0d+00
    di = 0.0d+00
    oneeff = 0.0d+00
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
        cmplxint = DCMPLX(one_elec_int_r(i, i), one_elec_int_i(i, i))
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

            i2r = inttwr(i, i, j, j)
            if (.not. realonly%is_realonly()) then
                i2i = inttwi(i, i, j, j)
            end if
            cmplxint = DCMPLX(i2r, i2i)
            energyHF(2) = energyHF(2) + cmplxint

            i2r = inttwr(i, j, j, i)
            if (.not. realonly%is_realonly()) then
                i2i = inttwi(i, j, j, i)
            end if
            cmplxint = DCMPLX(i2r, i2i)
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
        cmplxint = DCMPLX(one_elec_int_r(i, i), one_elec_int_i(i, i))
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
            if (.not. realonly%is_realonly()) then
                i2i = inttwi(i, i, j, j)
            end if
            cmplxint = DCMPLX(i2r, i2i)
            energy(iroot, 2) = energy(iroot, 2) + cmplxint

            i2r = inttwr(i, j, j, i)
            if (.not. realonly%is_realonly()) then
                i2i = inttwi(i, j, j, i)
            end if
            cmplxint = DCMPLX(i2r, i2i)
            energy(iroot, 2) = energy(iroot, 2) - cmplxint

        end do
    end do
    energy(iroot, 2) = 0.5d+00*energy(iroot, 2)
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
                if (.not. realonly%is_realonly()) then
                    i2i = inttwi(i, j, k, k)
                end if
                cmplxint = DCMPLX(i2r, i2i)
                oneeff = oneeff + cmplxint

                i2r = inttwr(i, k, k, j)
                if (.not. realonly%is_realonly()) then
                    i2i = inttwi(i, k, k, j)
                end if

                cmplxint = DCMPLX(i2r, i2i)
                oneeff = oneeff - cmplxint

            end do
            cmplxint = DCMPLX(one_elec_int_r(i, j), one_elec_int_i(i, j))
            oneeff = oneeff + cmplxint

            if (i == j) oneeff = oneeff*0.5d+00

            ii = convert_global_to_active_idx(i)
            jj = convert_global_to_active_idx(j)
            if (realonly%is_realonly()) then
                Call dim1_density_R(ii, jj, dr) ! dr is always zero
            else
                Call dim1_density(ii, jj, dr, di)
            end if
            dens = DCMPLX(dr, di)
            energy(iroot, 3) = energy(iroot, 3) + oneeff*dens

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

                    i2r = inttwr(i, j, k, l)
                    if (.not. realonly%is_realonly()) then
                        i2i = inttwi(i, j, k, l)
                    end if
                    cmplxint = DCMPLX(i2r, i2i)

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
                    energy(iroot, 4) = energy(iroot, 4) - dens*cmplxint

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

#ifdef HAVE_MPI
    call allreduce_wrapper(mat=energy(iroot, 1:4))
    call allreduce_wrapper(mat=energyHF(1:2))
#endif

    if (rank == 0) then
        print '(a,x,i0)', 'selectroot', iroot
        print *, 'core energy =', ecore
        print *, 'energy 1 =', energy(iroot, 1)
        print *, 'energy 2 =', energy(iroot, 2)
        print *, 'energy 3 =', energy(iroot, 3)
        print *, 'energy 4 =', energy(iroot, 4)
        print *, 't-energy(1-4)', sum(energy(iroot, :))
        print *, 't-energy', eigen(iroot) - ecore
        print *, 'C the error ', eigen(iroot) - ecore - sum(energy(iroot, :))
        print *, 'energy HF  =', sum(energyHF) + ecore
    end if
    Call memminus(KIND(energy), SIZE(energy), 1); deallocate (energy)
    if (rank == 0) print *, 'e0test end'
End subroutine e0test