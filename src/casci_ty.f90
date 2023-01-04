! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE casci_ty

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use module_file_manager
    use four_caspt2_module
    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer :: comb, j0, j, i0, irec, cimat
    real*8 :: thresd

    complex*16, allocatable :: mat(:, :)
    real*8, allocatable     :: ecas(:)
    logical                 :: cutoff
    character*20            :: filename
    real(8) :: expected_mem
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1

    ndet = comb(nact, nelec)
    if (rank == 0) print *, 'ndet', ndet
    Call casdet_ty
    if (rank == 0) then
        print *, "before allocate mat(ndet,ndet)"
        print '("Current Memory is ",F10.2,"MB")', tmem/1024/1024
        print *, 'kind of complex16 array named mat is ', kind(mat)
        expected_mem = tmem + (ndet**2)*16
        print *, 'expected used memory after allocate mat is ', expected_mem/1024/1024, 'MB'
    end if

    Allocate (mat(ndet, ndet)); Call memplus(KIND(mat), SIZE(mat), 2)
    if (rank == 0) print *, "end allocate mat(ndet,ndet)"
    Call casmat(mat)

    if (rank == 0) print *, 'before allocate ecas(ndet)'
    Allocate (ecas(ndet))
    if (rank == 0) print *, 'allocate ecas(ndet)'
    ecas = 0.0d+00
    thresd = 1.0d-15
    cutoff = .FALSE.
    if (rank == 0) print *, 'Start mat cdiag'
    datetmp1 = date0; datetmp0 = date0

    Call timing(date0, tsec0, datetmp0, tsectmp0)
    tsectmp1 = tsectmp0

    if (rank == 0) print *, 'ndet before cdiag', ndet
    Call cdiag(mat, ndet, ndet, ecas, thresd, cutoff)

    if (rank == 0) print *, 'End mat cdiag'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
! Print out CI matrix!
    if (rank == 0) then ! Only master ranks are allowed to create files used by CASPT2 except for MDCINTNEW.
        print *, 'debug1'
        cimat = default_unit
        filename = 'CIMAT'
        call open_unformatted_file(unit=cimat, file=filename, status='replace')
        write (cimat) ndet
        write (cimat) idet(1:ndet)
        write (cimat) ecas(1:ndet)
        write (cimat) 2**nact - 1 ! idetrの配列の要素数
        write (cimat) idetr(1:2**nact - 1)
        close (cimat)

! Print out C1 matrix!

! Print out CI matrix!

        print *, 'debug2'

        cimat = default_unit
        filename = 'CIMAT1'
        call open_unformatted_file(unit=cimat, file=filename, status='replace')
        write (cimat) ndet
        write (cimat) idet(1:ndet)
        write (cimat) ecas(1:ndet)
        write (cimat) mat(1:ndet, 1:ndet)
        close (cimat)
    end if
! Print out C1 matrix!

    if (rank == 0) print *, 'debug3'
    Allocate (cir(ndet, selectroot:selectroot)); Call memplus(KIND(cir), SIZE(cir), 1)
    Allocate (cii(ndet, selectroot:selectroot)); Call memplus(KIND(cii), SIZE(cii), 1)
    Allocate (eigen(nroot)); Call memplus(KIND(eigen), SIZE(eigen), 1)

    eigen(:) = 0.0d+00
    cir(:, :) = 0.0d+00
    cii(:, :) = 0.0d+00
    eigen(1:nroot) = ecas(1:nroot) + ecore
    cir(1:ndet, selectroot) = DBLE(mat(1:ndet, selectroot))
    cii(1:ndet, selectroot) = DIMAG(mat(1:ndet, selectroot))
    Deallocate (ecas)
    if (rank == 0) then
        print *, 'debug4'

        print '("CASCI ENERGY FOR ",I2," STATE")', totsym
        Do irec = 1, nroot
            print '(I4,F30.15)', irec, eigen(irec)
        End do
    end if
    do j = 1, ndet
        if (ABS(DIMAG(mat(j, selectroot))) > thres) then
            realcvec = .false.
        end if
    end do
    if (rank == 0) then
        do irec = 1, nroot
            print '("Root = ",I4)', irec
            do j = 1, ndet
                if ((ABS(mat(j, irec))**2) > 1.0d-02) then
                    i0 = idet(j)
                    print *, (btest(i0, j0), j0=0, nact - 1)
                    print '(I4,2(3X,E14.7)," Weights ",E14.7)', &
                    & j, mat(j, irec), &
                    & ABS(mat(j, irec))**2
                end if
            end do
        end do
    end if
    Deallocate (mat); Call memminus(KIND(mat), SIZE(mat), 2)
end subroutine casci_ty

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
