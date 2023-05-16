! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE fockcasci ! TO MAKE FOCK MATRIX for CASCI state

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables
    use module_index_utils, only: convert_global_to_active_idx
#ifdef HAVE_MPI
    use module_mpi
#endif
    Implicit NONE

    integer :: i, j, k, l
    integer :: kact, lact
    real(8) :: dr, di
    complex*16 :: dens
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!! NOW MAKE FOCK MATRIX FOR CASCI STATE
!! fij = hij + SIGUMA_kl[<0|Ekl|0>{(ij|kl)-(il|kj)}
    if (rank == 0) print *, 'enter building fock matrix'
    datetmp0 = initdate
    tsectmp0 = inittime
    call timing(datetmp0, tsectmp0, datetmp1, tsectmp1)
    datetmp0 = datetmp1
    tsectmp0 = tsectmp1

! Initialization
    dr = 0.0d+00; di = 0.0d+00; dens = 0.0d+00
    fock_cmplx(:, :) = 0.0d+00

    if (rank == 0) print *, 'enter building fock matrix'
!$OMP parallel private(i,j,k,l,dr,di,dens)
!$OMP do schedule(dynamic,2)
    do i = rank + 1, ninact + nact, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        do j = i, ninact + nact

            fock_cmplx(i, j) = DCMPLX(one_elec_int_r(i, j), one_elec_int_i(i, j))

            do k = 1, ninact
                fock_cmplx(i, j) = fock_cmplx(i, j) + DCMPLX(inttwr(i, j, k, k), inttwi(i, j, k, k))
                fock_cmplx(i, j) = fock_cmplx(i, j) - DCMPLX(inttwr(i, k, k, j), inttwi(i, k, k, j))
            End do
            do k = ninact + 1, ninact + nact              ! ACTIVE SPACE
                do l = ninact + 1, ninact + nact           ! ACTIVE SPACE
                    kact = convert_global_to_active_idx(k)
                    lact = convert_global_to_active_idx(l)
                    Call dim1_density(kact, lact, dr, di)
                    dens = DCMPLX(dr, di)
                    fock_cmplx(i, j) = fock_cmplx(i, j) + dens*DCMPLX(inttwr(i, j, k, l), inttwi(i, j, k, l))
                    fock_cmplx(i, j) = fock_cmplx(i, j) - dens*DCMPLX(inttwr(i, l, k, j), inttwi(i, l, k, j))
                End do
            End do

            fock_cmplx(j, i) = DCONJG(fock_cmplx(i, j))
        end do
    end do
!$OMP end do

!$OMP do schedule(dynamic,2)
    do i = ninact + nact + 1 + rank, ninact + nact + nsec, nprocs
        do j = i, ninact + nact + nsec
            fock_cmplx(i, j) = DCMPLX(one_elec_int_r(i, j), one_elec_int_i(i, j))
            do k = 1, ninact
                fock_cmplx(i, j) = fock_cmplx(i, j) + DCMPLX(int2r_f1(i, j, k, k), int2i_f1(i, j, k, k))
                fock_cmplx(i, j) = fock_cmplx(i, j) - DCMPLX(int2r_f2(i, k, k, j), int2i_f2(i, k, k, j))
            End do

            do k = ninact + 1, ninact + nact              ! ACTIVE SPACE
                do l = ninact + 1, ninact + nact           ! ACTIVE SPACE
                    kact = convert_global_to_active_idx(k)
                    lact = convert_global_to_active_idx(l)
                    Call dim1_density(kact, lact, dr, di)
                    dens = DCMPLX(dr, di)
                    fock_cmplx(i, j) = fock_cmplx(i, j) + dens*DCMPLX(int2r_f1(i, j, k, l), int2i_f1(i, j, k, l))
                    fock_cmplx(i, j) = fock_cmplx(i, j) - dens*DCMPLX(int2r_f2(i, l, k, j), int2i_f2(i, l, k, j))
                End do
            End do

            fock_cmplx(j, i) = DCONJG(fock_cmplx(i, j))
        end do
    end do
!$OMP end do
!$OMP end parallel
    if (rank == 0) print *, 'fockcasci before fock_cmplx allreduce'
    call timing(datetmp0, tsectmp0, datetmp1, tsectmp1)
    datetmp0 = datetmp1
    tsectmp0 = tsectmp1
#ifdef HAVE_MPI
    call allreduce_wrapper(mat=fock_cmplx(1:nmo, 1:nmo))
#endif
    if (rank == 0) print *, 'fockcasci end'
end subroutine fockcasci

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE fockcasci_real ! TO MAKE FOCK MATRIX for CASCI state

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables
    use module_index_utils, only: convert_global_to_active_idx
#ifdef HAVE_MPI
    use module_mpi
#endif
    Implicit NONE

    integer :: j, i, k, l
    integer :: kact, lact
    real(8) :: dr
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!! NOW MAKE FOCK MATRIX FOR CASCI STATE
!! fij = hij + SIGUMA_kl[<0|Ekl|0>{(ij|kl)-(il|kj)}
    if (rank == 0) print *, 'enter building fock matrix'
    datetmp0 = initdate
    tsectmp0 = inittime
    call timing(datetmp0, tsectmp0, datetmp1, tsectmp1)
    datetmp0 = datetmp1
    tsectmp0 = tsectmp1

! Initialization
    dr = 0.0d+00
    fock_real(:, :) = 0.0d+00

    if (rank == 0) print *, 'enter building fock matrix'
!$OMP parallel private(i,j,k,l,dr)
!$OMP do schedule(dynamic,2)
    do i = rank + 1, ninact + nact, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        do j = i, ninact + nact

            fock_real(i, j) = one_elec_int_r(i, j)
            do k = 1, ninact
                fock_real(i, j) = fock_real(i, j) + inttwr(i, j, k, k)
                fock_real(i, j) = fock_real(i, j) - inttwr(i, k, k, j)
            End do
            do k = ninact + 1, ninact + nact              ! ACTIVE SPACE
                do l = ninact + 1, ninact + nact           ! ACTIVE SPACE

                    kact = convert_global_to_active_idx(k)
                    lact = convert_global_to_active_idx(l)
                    Call dim1_density_R(kact, lact, dr)
                    fock_real(i, j) = fock_real(i, j) + dr*inttwr(i, j, k, l)
                    fock_real(i, j) = fock_real(i, j) - dr*inttwr(i, l, k, j)

                End do
            End do

            fock_real(j, i) = fock_real(i, j)
        end do
    end do
!$OMP end do

!$OMP do schedule(dynamic,2)
    do i = ninact + nact + 1 + rank, ninact + nact + nsec, nprocs
        do j = i, ninact + nact + nsec
            fock_real(i, j) = one_elec_int_r(i, j)

            do k = 1, ninact
                fock_real(i, j) = fock_real(i, j) + int2r_f1(i, j, k, k)
                fock_real(i, j) = fock_real(i, j) - int2r_f2(i, k, k, j)
            End do
            do k = ninact + 1, ninact + nact              ! ACTIVE SPACE
                do l = ninact + 1, ninact + nact           ! ACTIVE SPACE
                    kact = convert_global_to_active_idx(k)
                    lact = convert_global_to_active_idx(l)
                    Call dim1_density_R(kact, lact, dr)
                    fock_real(i, j) = fock_real(i, j) + dr*int2r_f1(i, j, k, l)
                    fock_real(i, j) = fock_real(i, j) - dr*int2r_f2(i, l, k, j)

                End do
            End do
            fock_real(j, i) = fock_real(i, j)
        end do
    end do
!$OMP end do
!$OMP end parallel
    if (rank == 0) print *, 'fockcasci_real before fock_real allreduce'
    call timing(datetmp0, tsectmp0, datetmp1, tsectmp1)
    datetmp0 = datetmp1
    tsectmp0 = tsectmp1
#ifdef HAVE_MPI
    call allreduce_wrapper(mat=fock_real(1:nmo, 1:nmo))
#endif
    if (rank == 0) print *, 'fockcasci_real end'
end subroutine fockcasci_real
