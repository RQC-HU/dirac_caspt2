! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE fockhf1_ty ! TO CALCULATE FOCK MATRIX OF HF STATE, A TEST

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer :: j, i, k, n

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    debug = .TRUE.
    thres = 1.0d-15
!        thres = 0.0d+00

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

!! TEST TO CALCULATE FOCK MATRIX OF HF STATE fpq = hpq + SIGUMA_r[(pq|rr)-(pr|qr)]
!! THIS MUST BE DIAGONAL MATRIX AND DIAGONAL ELEMENTS CORESPONDS TO SPINOR ENERGIES.
    if (rank == 0) then
        print *, ' '
        print *, 'FOR TEST, FOCK MATRIX OF HF STATE IS CALCULATED '
    end if
    n = 0
    f = 0.0d+00

    !$OMP parallel do private(j,k)
    do i = rank + 1, ninact + nact, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        do j = i, ninact + nact
            f(i, j) = DCMPLX(oner(i, j), onei(i, j))
            do k = 1, ninact + nelec

                f(i, j) = f(i, j) + CMPLX(inttwr(i, j, k, k), inttwi(i, j, k, k), 16)
                f(i, j) = f(i, j) - CMPLX(inttwr(i, k, k, j), inttwi(i, k, k, j), 16)

            End do           ! k

            f(j, i) = DCONJG(f(i, j))
        End do       ! j
    End do          ! i

    !$OMP parallel do private(j,k)
    do i = rank + ninact + nact + 1, ninact + nact + nsec, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        do j = i, ninact + nact + nsec
            f(i, j) = DCMPLX(oner(i, j), onei(i, j))
            do k = 1, ninact + nelec

                f(i, j) = f(i, j) + DCMPLX(int2r_f1(i, j, k, k), int2i_f1(i, j, k, k))
                f(i, j) = f(i, j) - DCMPLX(int2r_f2(i, k, k, j), int2i_f2(i, k, k, j))

            End do           ! k

            f(j, i) = DCONJG(f(i, j))

        End do       ! j
    End do          ! i
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, f(1, 1), nmo**2, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    if (rank == 0) then
        print *, ' '
        print *, 'OFF DIAGONAL ELEMENTS OF FOCK MATRIX WHICH IS LARGER THAN 1.0d-06 '
        print *, ' '
        do i = 1, ninact + nact + nsec
            do j = i, ninact + nact + nsec
                if ((i /= j) .and. (ABS(f(i, j)) > 1.0d-6)) then
                    print '(2I4,2E20.10)', i, j, f(i, j)
                end if
            end do
        end do
        print *, ' '
        print *, 'THESE DIAGONAL ELEMENTS SHOULD BE CORESPOND TO HF SPINOR ENERGY '
        print *, ' '
        print *, '  NO.   Spinor Energy(Re)   Spinor Energy(Im) '&
        &, 'Spinor Energy (HF)        ERROR'
        do i = 1, ninact + nact + nsec
            print '(I4,4E20.10)', i, f(i, i), orbmo(i), orbmo(i) - dble(f(i, i))
        end do

        print *, 'fockhf end'
    end if
end SUBROUTINE fockhf1_ty
