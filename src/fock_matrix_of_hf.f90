! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE fock_matrix_of_hf_complex ! TO CALCULATE FOCK MATRIX OF HF STATE, A TEST

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables
#ifdef HAVE_MPI
    use module_mpi
#endif
    Implicit NONE

    integer :: j, i, k, n

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    debug = .TRUE.

!! TEST TO CALCULATE FOCK MATRIX OF HF STATE fpq = hpq + SIGUMA_r[(pq|rr)-(pr|qr)]
!! THIS MUST BE DIAGONAL MATRIX AND DIAGONAL ELEMENTS CORESPONDS TO SPINOR ENERGIES.
    if (rank == 0) then
        print *, ' '
        print *, 'FOR TEST, FOCK MATRIX OF HF STATE IS CALCULATED '
    end if
    n = 0
    fock_cmplx = 0.0d+00

!$OMP parallel private(i,j,k)
!$OMP do
    do i = rank + 1, global_act_end, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        do j = i, global_act_end
            fock_cmplx(i, j) = DCMPLX(one_elec_int_r(i, j), one_elec_int_i(i, j))
            do k = 1, ninact + nelec

                fock_cmplx(i, j) = fock_cmplx(i, j) + DCMPLX(inttwr(i, j, k, k), inttwi(i, j, k, k))
                fock_cmplx(i, j) = fock_cmplx(i, j) - DCMPLX(inttwr(i, k, k, j), inttwi(i, k, k, j))

            End do           ! k

            fock_cmplx(j, i) = DCONJG(fock_cmplx(i, j))
        End do       ! j
    End do          ! i
!$OMP end do
!$OMP do
    do i = rank + global_sec_start, global_sec_end, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        do j = i, global_sec_end
            fock_cmplx(i, j) = DCMPLX(one_elec_int_r(i, j), one_elec_int_i(i, j))
            do k = 1, ninact + nelec

                fock_cmplx(i, j) = fock_cmplx(i, j) + DCMPLX(int2r_f1(i, j, k, k), int2i_f1(i, j, k, k))
                fock_cmplx(i, j) = fock_cmplx(i, j) - DCMPLX(int2r_f2(i, k, k, j), int2i_f2(i, k, k, j))

            End do           ! k

            fock_cmplx(j, i) = DCONJG(fock_cmplx(i, j))

        End do       ! j
    End do          ! i
!$OMP end do
!$OMP end parallel

#ifdef HAVE_MPI
    call allreduce_wrapper(mat=fock_cmplx(1:nmo, 1:nmo))
#endif

    if (rank == 0) then
        print *, ' '
        print *, 'OFF DIAGONAL ELEMENTS OF FOCK MATRIX WHICH IS LARGER THAN 1.0d-06 '
        print *, ' '
        do i = 1, global_sec_end
            do j = i, global_sec_end
                if ((i /= j) .and. (ABS(fock_cmplx(i, j)) > 1.0d-6)) then
                    print '(2I4,2E20.10)', i, j, fock_cmplx(i, j)
                end if
            end do
        end do
        print *, ' '
        print *, 'THESE DIAGONAL ELEMENTS SHOULD BE CORESPOND TO HF SPINOR ENERGY '
        print *, ' '
        print *, '  NO.   Spinor Energy(Re)   Spinor Energy(Im) '&
        &, 'Spinor Energy (HF)        ERROR'
        do i = 1, global_sec_end
            print '(I4,4E20.10)', i, fock_cmplx(i, i), caspt2_mo_energy(i), caspt2_mo_energy(i) - dble(fock_cmplx(i, i))
        end do

        print *, 'fockhf end'
    end if
end SUBROUTINE fock_matrix_of_hf_complex
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE fock_matrix_of_hf_real ! TO CALCULATE FOCK MATRIX OF HF STATE, A TEST

    ! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    ! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables
#ifdef HAVE_MPI
    use module_mpi
#endif
    Implicit NONE

    integer :: j, i, k, n

    ! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    ! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    debug = .TRUE.

    !! TEST TO CALCULATE FOCK MATRIX OF HF STATE fpq = hpq + SIGUMA_r[(pq|rr)-(pr|qr)]
    !! THIS MUST BE DIAGONAL MATRIX AND DIAGONAL ELEMENTS CORESPONDS TO SPINOR ENERGIES.
    if (rank == 0) then
        print *, ' '
        print *, 'FOR TEST, FOCK MATRIX OF HF STATE IS CALCULATED '
    end if
    n = 0
    fock_real = 0.0d+00

!$OMP parallel private(i,j,k)
!$OMP do
    do i = rank + 1, global_act_end, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        do j = i, global_act_end
            fock_real(i, j) = one_elec_int_r(i, j)
            do k = 1, ninact + nelec

                fock_real(i, j) = fock_real(i, j) + inttwr(i, j, k, k)
                fock_real(i, j) = fock_real(i, j) - inttwr(i, k, k, j)

            End do           ! k

            fock_real(j, i) = fock_real(i, j)
        End do       ! j
    End do          ! i
!$OMP end do
!$OMP do
    do i = rank + global_sec_start, global_sec_end, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        do j = i, global_sec_end
            fock_real(i, j) = one_elec_int_r(i, j)
            do k = 1, ninact + nelec

                fock_real(i, j) = fock_real(i, j) + int2r_f1(i, j, k, k)
                fock_real(i, j) = fock_real(i, j) - int2r_f2(i, k, k, j)

            End do           ! k

            fock_real(j, i) = fock_real(i, j)

        End do       ! j
    End do          ! i
!$OMP end do
!$OMP end parallel

#ifdef HAVE_MPI
    call allreduce_wrapper(mat=fock_real(1:nmo, 1:nmo))
#endif

    if (rank == 0) then
        print *, ' '
        print *, 'OFF DIAGONAL ELEMENTS OF FOCK MATRIX WHICH IS LARGER THAN 1.0d-06 '
        print *, ' '
        do i = 1, global_sec_end
            do j = i, global_sec_end
                if ((i /= j) .and. (ABS(fock_real(i, j)) > 1.0d-6)) then
                    print '(2I4,2E20.10)', i, j, fock_real(i, j)
                end if
            end do
        end do
        print *, ' '
        print *, 'THESE DIAGONAL ELEMENTS SHOULD BE CORESPOND TO HF SPINOR ENERGY '
        print *, ' '
        print *, '  NO.   Spinor Energy(Re)   Spinor Energy(Im) '&
        &, 'Spinor Energy (HF)        ERROR'
        do i = 1, global_sec_end
            print '(I4,4E20.10)', i, fock_real(i, i), caspt2_mo_energy(i), caspt2_mo_energy(i) - dble(fock_real(i, i))
        end do

        print *, 'fockhf end'
    end if
end SUBROUTINE fock_matrix_of_hf_real
