! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE fockcasci_ty ! TO MAKE FOCK MATRIX for CASCI state

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module
#ifdef HAVE_MPI
    use module_mpi
#endif
    Implicit NONE

    integer :: j, i, k, l
    real*8 :: dr, di
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

    f = 0.0d+00

    if (rank == 0) print *, 'enter building fock matrix'
!$OMP parallel private(i,j,k,l,dr,di,dens)
!$OMP do schedule(dynamic,2)
    do i = rank + 1, ninact + nact, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        do j = i, ninact + nact

            f(i, j) = DCMPLX(one_elec_int_r(i, j), one_elec_int_i(i, j))
            do k = 1, ninact
                f(i, j) = f(i, j) + DCMPLX(inttwr(i, j, k, k), inttwi(i, j, k, k))
                f(i, j) = f(i, j) - DCMPLX(inttwr(i, k, k, j), inttwi(i, k, k, j))
            End do           ! k
            do k = ninact + 1, ninact + nact              ! ACTIVE SPACE
                do l = ninact + 1, ninact + nact           ! ACTIVE SPACE
                    If (realcvec) then
                        Call dim1_density_R(k - ninact, l - ninact, dr)
                        f(i, j) = f(i, j) + dr*DCMPLX(inttwr(i, j, k, l), inttwi(i, j, k, l))
                        f(i, j) = f(i, j) - dr*DCMPLX(inttwr(i, l, k, j), inttwi(i, l, k, j))
                    Else
                        dr = 0.0d+00
                        Call dim1_density(k - ninact, l - ninact, dr, di)
                        dens = DCMPLX(dr, di)

                        f(i, j) = f(i, j) + dens*DCMPLX(inttwr(i, j, k, l), inttwi(i, j, k, l))
                        f(i, j) = f(i, j) - dens*DCMPLX(inttwr(i, l, k, j), inttwi(i, l, k, j))
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
            f(i, j) = DCMPLX(one_elec_int_r(i, j), one_elec_int_i(i, j))
            do k = 1, ninact
                f(i, j) = f(i, j) + DCMPLX(int2r_f1(i, j, k, k), int2i_f1(i, j, k, k))
                f(i, j) = f(i, j) - DCMPLX(int2r_f2(i, k, k, j), int2i_f2(i, k, k, j))
            End do           ! k

            do k = ninact + 1, ninact + nact              ! ACTIVE SPACE
                do l = ninact + 1, ninact + nact           ! ACTIVE SPACE
                    If (realcvec) then
                        Call dim1_density_R(k - ninact, l - ninact, dr)
                        f(i, j) = f(i, j) + dr*DCMPLX(int2r_f1(i, j, k, l), int2i_f1(i, j, k, l))
                        f(i, j) = f(i, j) - dr*DCMPLX(int2r_f2(i, l, k, j), int2i_f2(i, l, k, j))
                    Else
                        Call dim1_density(k - ninact, l - ninact, dr, di)
                        dens = DCMPLX(dr, di)
                        f(i, j) = f(i, j) + dens*DCMPLX(int2r_f1(i, j, k, l), int2i_f1(i, j, k, l))
                        f(i, j) = f(i, j) - dens*DCMPLX(int2r_f2(i, l, k, j), int2i_f2(i, l, k, j))
                    End if
                End do     ! l
            End do        ! k

            f(j, i) = DCONJG(f(i, j))
        end do       ! j
    end do          ! i
!$OMP end do
!$OMP end parallel
    if (rank == 0) print *, 'fockcasci before f allreduce'
    call timing(datetmp0, tsectmp0, datetmp1, tsectmp1)
    datetmp0 = datetmp1
    tsectmp0 = tsectmp1
#ifdef HAVE_MPI
    call allreduce_wrapper(mat=f(1:nmo, 1:nmo))
#endif
    if (rank == 0) print *, 'fockcasci end'
end subroutine fockcasci_ty
