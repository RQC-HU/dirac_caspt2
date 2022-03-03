! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE casci_ty(totsym)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module
    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer, intent(in) :: totsym
    integer :: nbitsa, comb
    integer :: j0, j, i, i0, i1
    integer :: k0, l0, ii, jj, kk, ll, irec, cimat
    real*8 :: thresd

    complex*16, allocatable :: mat(:, :)
    real*8, allocatable     :: ecas(:)
    logical                 :: cutoff
    character*20            :: filename
    real(8) :: expected_mem
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1
#ifdef CASMAT_DEBUG
    character(50)           :: chr_rank, matfilename
    integer, parameter      :: mat_unit_num = 1000
    real*8             :: matvalr, matvali
    complex(16), allocatable :: mat_not_zero(:)
    integer, allocatable    :: mat_i_order(:), mat_j_order(:)
    integer                 :: not_zero_sum, not_zero_count, loop_idx, read_count, loop_sum, not_zero_count_modified
    integer :: idx1, idx2, matcount, matdiagcount
    real(8) :: not_zero_percentage
#endif

    ndet = comb(nact, nelec)
    if (rank == 0) then ! Process limits for output
        write (normaloutput, *) 'ndet', ndet
    end if
    Call casdet_ty(totsym)
    if (rank == 0) then
        write (normaloutput, *) "before allocate mat(ndet,ndet)"
        write (normaloutput, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024
        write (normaloutput, *) 'kind of complex16 array named mat is ', kind(mat)
        expected_mem = tmem + (ndet**2)*16
        write (normaloutput, *) 'expected used memory after allocate mat is ', expected_mem/1024/1024, 'MB'
    end if

#ifdef CASMAT_DEBUG
    not_zero_count = 0
    call casmat_modified(not_zero_count)
    if (rank == 0) write (normaloutput, *) 'Noda end casmat_modified'
    not_zero_sum = not_zero_count
    call MPI_Allreduce(MPI_IN_PLACE, not_zero_sum, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)

    if (rank == 0) then
        not_zero_percentage = (real(not_zero_sum, 8)/real((ndet**2), 8))*100
        write (normaloutput, *) 'Noda end casmat_modified, not_zero_sum:', not_zero_sum, &
            ' all', ndet**2, ' percentage', not_zero_percentage, '%'
    end if
    if (rank == 0) then
        allocate (mat_not_zero(not_zero_sum))
        allocate (mat_i_order(not_zero_sum))
        allocate (mat_j_order(not_zero_sum))
        if (rank == 0) write (normaloutput, *) "Noda nprocs check", nprocs
        read_count = 0
        not_zero_count_modified = 0
        do loop_idx = 0, nprocs - 1
            write (chr_rank, *) loop_idx
            matfilename = 'mat'//trim(adjustl(chr_rank))
            if (rank == 0) write (normaloutput, *) "Noda matfilename check", trim(matfilename)
            open (mat_unit_num, file=matfilename, form='unformatted')
103         read (mat_unit_num, err=101, end=102) i, j, matvalr, matvali
            if (i == j) then
                read_count = read_count + 1
                not_zero_count_modified = not_zero_count_modified + 1
            else
                read_count = read_count + 1
                not_zero_count_modified = not_zero_count_modified + 2
            end if
            mat_i_order(read_count) = i
            mat_j_order(read_count) = j
            mat_not_zero(read_count) = matvalr

            goto 103
101         write (*, *) 'mat_read_error', rank, nprocs
102         if (rank == 0) then
                not_zero_percentage = (real(not_zero_count_modified, 8)/real((ndet**2), 8))*100
                write (normaloutput, *) 'Noda end read casmat_modified, not_zero_count_modified:', not_zero_count_modified, &
                    ' all', ndet**2, ' percentage', not_zero_percentage, '%'
            end if
            close (mat_unit_num)
        end do
    end if
    loop_sum = read_count
    ! if (rank == 0) then
    !     write (normaloutput, *) 'mat mat_not_zero diff check'
    !     do loop_idx = 1, not_zero_sum
    !         if (mat(mat_i_order(loop_idx), mat_j_order(loop_idx)) /= mat_not_zero(loop_idx)) then
    !             write (normaloutput, '(A,I5,A,I5,A,E20.10)') 'Noda mat/=mat_not_zero i: ', mat_i_order(loop_idx), &
    !                 ' j: ', mat_j_order(loop_idx), 'mat ', mat(mat_i_order(loop_idx), mat_j_order(loop_idx)), mat_not_zero(loop_idx)
    !         end if
    !         write (normaloutput, '(A,I5,A,I5,A,2E20.10)') 'Noda modified mat i: ', mat_i_order(loop_idx), ' j: ',&
    !                 mat_j_order(loop_idx), ' mat:', mat_not_zero(loop_idx)
    !     end do
    !     write (normaloutput, *) 'end mat mat_not_zero diff check'
    ! end if
#endif
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
#ifdef BIG_MAT
    if (rank == 0) then
        Allocate (mat(ndet, ndet)); Call memplus(KIND(mat), SIZE(mat), 2)
        write (normaloutput, *) "end allocate mat(ndet,ndet)"
        write (normaloutput, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024
        Call casmat(mat)
    end if
#else
    Allocate (mat(ndet, ndet)); Call memplus(KIND(mat), SIZE(mat), 2)
    if (rank == 0) then
        write (normaloutput, *) "end allocate mat(ndet,ndet)"
        write (normaloutput, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024
    end if
    Call casmat(mat)
#endif
#ifdef CASMAT_DEBUG
    if (rank == 0) then
        matcount = 0
        matdiagcount = 0
        write (normaloutput, *) 'Noda end casmat'
        do idx1 = 1, ndet
            do idx2 = idx1, ndet
                if (mat(idx2, idx1) /= 0.0d+00) then
                    if (idx1 == idx2) then
                        matdiagcount = matdiagcount + 1
                        matcount = matcount + 1
                    else
                        matcount = matcount + 2
                    end if
                    ! write (normaloutput, *) idx2, idx1, mat(idx2, idx1), mat(idx1, idx2)
                end if
            end do
        end do
        not_zero_percentage = (real(matcount, 8)/real((ndet**2), 8))*100
        write (normaloutput, *) 'matcount', matcount, ' all', ndet**2, ' percentage', not_zero_percentage, '%'
        write (normaloutput, *) 'matdiagcount', matdiagcount, ' ndet', ndet
    end if
    if (rank == 0) write (normaloutput, *) 'mat,mat_not_zero diff check'
    do read_count = 1, loop_sum
!   if (rank == 0.and.abs(real(mat(mat_i_order(read_count), mat_j_order(read_count)),8)-real(mat_not_zero(read_count),8))>1.0d-15)then
       if (rank == 0 .and. real(mat(mat_i_order(read_count), mat_j_order(read_count)), 8) /= real(mat_not_zero(read_count), 8)) then
            write (normaloutput, '(A,I5,A,I5,A,2E20.10,A,E20.10)') 'Noda check i: ', mat_i_order(read_count), &
                ' j: ', mat_j_order(read_count), &
                'mat ', real(mat(mat_i_order(read_count), mat_j_order(read_count)), 8), real(mat_not_zero(read_count), 8), &
                'diff ', (real(mat(mat_i_order(read_count), mat_j_order(read_count)), 8) - real(mat_not_zero(read_count), 8))
        end if
    end do
    if (rank == 0) write (normaloutput, *) 'mat,mat_not_zero diff check end'
#endif
    if (rank == 0) then ! Process limits for output
        write (normaloutput, *) 'before allocate ecas(ndet)'
    end if
    Allocate (ecas(ndet))
    if (rank == 0) then ! Process limits for output
        write (normaloutput, *) 'allocate ecas(ndet)'
    end if
    ecas = 0.0d+00
    thresd = 1.0d-15
    cutoff = .FALSE.
    if (rank == 0) write (normaloutput, *) 'Start mat cdiag'
    datetmp1 = date0; datetmp0 = date0

    Call timing(date0, tsec0, datetmp0, tsectmp0)
    tsectmp1 = tsectmp0

    if (rank == 0) write (normaloutput, *) 'Noda ndet before cdiag', ndet
#ifdef BIG_MAT
    if (rank == 0) then
        ! Only the master process has a matrix named mat.
        ! So only the master process will diagonalize this matrix.
        Call cdiag(mat, ndet, ndet, ecas, thresd, cutoff)
    end if
#else
    Call cdiag(mat, ndet, ndet, ecas, thresd, cutoff)
#endif
#ifdef CASMAT_DEBUG
    if (rank == 0) then
        write (normaloutput, *) 'Noda ndet after cdiag', ndet
        write (normaloutput, *) 'Noda cidag mat diag values'
        do idx1 = 1, ndet
            ! if (real(mat(idx1, idx1)) > 1.0d-15) then
            write (normaloutput, *) idx1, mat(idx1, idx1)
            ! end if
        end do
    end if
#endif
!    if (rank == 0) then
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) write (normaloutput, *) 'End mat cdiag'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
! Print out CI matrix!
    if (rank == 0) then ! Only master ranks are allowed to create files used by CASPT2 except for MDCINTNEW.
        write (normaloutput, *) 'debug1'
        cimat = 10
        filename = 'CIMAT'
        open (10, file='CIMAT', status='unknown', form='unformatted')
        write (10) ndet
        write (10) idet(1:ndet)
        write (10) ecas(1:ndet)
        write (10) 2**nact - 1 ! idetrの配列の要素数
        write (10) idetr(1:2**nact - 1)
!        write(10) mat(1:ndet,1:ndet)
        close (10)

! Print out C1 matrix!

! Print out CI matrix!

        write (normaloutput, *) 'debug2'

        cimat = 10
        filename = 'CIMAT1'
        open (10, file='CIMAT1', status='unknown', form='unformatted')
        write (10) ndet
        write (10) idet(1:ndet)
        write (10) ecas(1:ndet)
        write (10) mat(1:ndet, 1:ndet)
        close (10)
    end if
! Print out C1 matrix!

    if (rank == 0) then ! Process limits for output
        write (normaloutput, *) 'debug3'
    end if
    Allocate (cir(ndet, selectroot:selectroot)); Call memplus(KIND(cir), SIZE(cir), 1)
    Allocate (cii(ndet, selectroot:selectroot)); Call memplus(KIND(cii), SIZE(cii), 1)
    Allocate (eigen(nroot)); Call memplus(KIND(eigen), SIZE(eigen), 1)

    eigen(:) = 0.0d+00
    cir(:, :) = 0.0d+00
    cii(:, :) = 0.0d+00
#ifdef BIG_MAT
    if (rank == 0) then
        eigen(1:nroot) = ecas(1:nroot) + ecore
        cir(1:ndet, selectroot) = DBLE(mat(1:ndet, selectroot))
        cii(1:ndet, selectroot) = DIMAG(mat(1:ndet, selectroot))
    end if
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Bcast(ndet, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
    write (*, *) 'ndet,rank', ndet, rank
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    write (*, *) 'selectroot,rank', selectroot, rank
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    write (*, *) 'nroot,rank', nroot, rank
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Bcast(eigen(1), nroot, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    write (*, *) 'eigen rank', rank, eigen
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Bcast(cir(1, selectroot), ndet, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    write (*, *) 'cir rank', rank, cir(1, selectroot)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Bcast(cii(1, selectroot), ndet, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    write (*, *) 'cii rank', rank, cir(1, selectroot)
#endif
#else
    eigen(1:nroot) = ecas(1:nroot) + ecore
    cir(1:ndet, selectroot) = DBLE(mat(1:ndet, selectroot))
    cii(1:ndet, selectroot) = DIMAG(mat(1:ndet, selectroot))
#endif
#ifdef CASMAT_DEBUG
    do loop_idx = 1, ndet
        write (normaloutput, *) 'Noda cir check', loop_idx, real(mat(loop_idx, selectroot), 8)
    end do
#endif
    Deallocate (ecas)
    if (rank == 0) then ! Process limits for output
        write (normaloutput, *) 'debug4'

        write (normaloutput, '("CASCI ENERGY FOR ",I2," STATE")') totsym
        Do irec = 1, nroot
            write (normaloutput, '(I4,F30.15)') irec, eigen(irec)
        End do

        do j = 1, ndet
            if (ABS(DIMAG(mat(j, selectroot))) > thres) then
                realcvec = .false.
            end if
        end do

        do irec = 1, nroot
            write (normaloutput, '("Root = ",I4)') irec
            do j = 1, ndet
                if ((ABS(mat(j, irec))**2) > 1.0d-02) then
                    i0 = idet(j)
                    write (normaloutput, *) (btest(i0, j0), j0=0, nact - 1)
                    write (normaloutput, '(I4,2(3X,E14.7)," Weights ",E14.7)') &
                    & j, mat(j, irec), &
                    & ABS(mat(j, irec))**2
                end if
            end do
        end do
    end if
#ifdef BIG_MAT
    ! Only the master process has a matrix named mat.
    ! So the only master process deallocate this.
    if (rank == 0) then
        Deallocate (mat); Call memminus(KIND(mat), SIZE(mat), 2)
    end if
#else
    Deallocate (mat); Call memminus(KIND(mat), SIZE(mat), 2)
#endif
1000 end subroutine casci_ty

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

FUNCTION comb(n, m) RESULT(res)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    Implicit NONE

    integer :: n, m, i, j, res, m0

    j = 1

    if (n - m < m) then
        m0 = n - m
    else
        m0 = m
    end if

    Do i = n - m0 + 1, n
        j = j*i
    End do

    Do i = 1, m0
        j = j/i
    End do

    res = j
1000 end function comb
