! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE casci

! CASCI calculation

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use module_file_manager
    use four_caspt2_module
    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer :: comb, j0, j, i0, irec, unit_cimat
    real*8 :: cutoff_threshold

    complex*16, allocatable :: mat(:, :)
    real*8, allocatable     :: ecas(:)
    character*20            :: filename, chr_root
    real(8) :: expected_mem
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1

    ndet = comb(nact, nelec)
    if (rank == 0) print *, 'ndet', ndet
    Call casdet
    if (rank == 0) then
        print *, "before allocate mat(ndet,ndet)"
        print '("Current Memory is ",F10.2,"MB")', tmem/1024/1024
        print *, 'kind of complex16 array named mat is ', kind(mat)
        expected_mem = tmem + (ndet**2)*16
        print *, 'expected used memory after allocate mat is ', expected_mem/1024/1024, 'MB'
    end if

    ! Create a matrix for CI
    Allocate (mat(ndet, ndet)); Call memplus(KIND(mat), SIZE(mat), 2)
    if (rank == 0) print *, "end allocate mat(ndet,ndet)"
    Call casmat(mat)
    Allocate (ecas(ndet))
    ecas = 0.0d+00
    datetmp1 = date0; datetmp0 = date0
    Call timing(date0, tsec0, datetmp0, tsectmp0)
    tsectmp1 = tsectmp0

    ! Diagonalize the CI matrix
    if (rank == 0) then
        print *, 'Start mat cdiag'
        print *, 'ndet before cdiag', ndet
    end if
    cutoff_threshold = 0  ! No need to resolve linear dependence
    Call cdiag(mat, ndet, ndet, ecas, cutoff_threshold)
    if (rank == 0) print *, 'End mat cdiag'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0

! Print out CI matrix!
    if (rank == 0) then ! Only master ranks are allowed to create files used by CASPT2 except for MDCINTNEW.
        filename = 'CIMAT'
        call open_unformatted_file(unit=unit_cimat, file=filename, status='replace')
        write (unit_cimat) ndet
        write (unit_cimat) cas_idx(1:ndet)
        write (unit_cimat) ecas(1:ndet)
        write (unit_cimat) 2**nact - 1 ! cas_idx_reverseの配列の要素数
        write (unit_cimat) cas_idx_reverse(1:2**nact - 1)
        close (unit_cimat)

        filename = 'CIMAT1'
        call open_unformatted_file(unit=unit_cimat, file=filename, status='replace')
        write (unit_cimat) ndet
        write (unit_cimat) cas_idx(1:ndet)
        write (unit_cimat) ecas(1:ndet)
        write (unit_cimat) mat(1:ndet, 1:ndet)
        close (unit_cimat)
    end if

    Allocate (cir(ndet, selectroot:selectroot)); Call memplus(KIND(cir), SIZE(cir), 1)
    Allocate (cii(ndet, selectroot:selectroot)); Call memplus(KIND(cii), SIZE(cii), 1)
    Allocate (eigen(nroot)); Call memplus(KIND(eigen), SIZE(eigen), 1)

    ! Print out the results
    eigen(:) = 0.0d+00
    cir(:, :) = 0.0d+00
    cii(:, :) = 0.0d+00
    eigen(1:nroot) = ecas(1:nroot) + ecore
    cir(1:ndet, selectroot) = DBLE(mat(1:ndet, selectroot))
    cii(1:ndet, selectroot) = DIMAG(mat(1:ndet, selectroot))
    Deallocate (ecas)
    if (rank == 0) then
        print '("CASCI ENERGY FOR ",I2," STATE")', totsym
        Do irec = 1, nroot
            write (chr_root, '(I4)') irec
            print '("CASCI Total Energy ROOT",a,F30.15," a.u.")', trim(adjustl(chr_root)), eigen(irec)
        End do
    end if
    do j = 1, ndet
        if (ABS(DIMAG(mat(j, selectroot))) > global_threshold) then
            realcvec = .false.
        end if
    end do
    if (rank == 0) then
        do irec = 1, nroot
            print '("Root = ",I4)', irec
            do j = 1, ndet
                if ((ABS(mat(j, irec))**2) > 1.0d-02) then
                    i0 = cas_idx(j)
                    print *, (btest(i0, j0), j0=0, nact - 1)
                    print '(I4,2(3X,E14.7)," Weights ",E14.7)', &
                    & j, mat(j, irec), &
                    & ABS(mat(j, irec))**2
                end if
            end do
        end do
    end if
    Deallocate (mat); Call memminus(KIND(mat), SIZE(mat), 2)
end subroutine casci

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
end function comb
