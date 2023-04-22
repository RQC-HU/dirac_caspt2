! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE fockivo_co ! TO MAKE FOCK MATRIX for IVO

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module
    use module_file_manager

    Implicit NONE

    integer      :: j, i, k, i0, j0
    integer      :: isym, nv, numh
    integer      :: imo, iao, IMAX, iostat
    integer      :: unit_dfpcmo, unit_itrfmog, unit_itrfmou, unit_buf
    real*8       :: thresd
    logical      :: cutoff
    complex*16, allocatable  :: fsym(:, :) ! Symmetrized fock_ivo_matrix for particular irrep
    complex*16, allocatable  :: coeff(:, :)
    real*8, allocatable      :: BUF(:)  ! One dimensional array representing MO coeff. read from DFPCMO
    real*8, allocatable      :: wsym(:), eval(:)
    integer, allocatable     :: mosym(:)
    character*150 :: line0, line1, line2, line3, line4, line5

    ! for new code of IVO
    integer :: npg, neg, nbasg ! number of positronic gerade(g), electronic g, basis set for g
    integer :: npu, neu, nbasu ! number of positronic ungerade(u), electronic u, basis set for u
    ! integer :: nvcutg, nvcutu ! number of virtual cut (g) and (u).
    integer :: nsum !nsum = npg + neg + npu + neu
    integer :: nv0, ngu, A, B ! A and B are dammy indices written in DFPCMO, A is nfsym in DIRAC
    integer, allocatable :: syminfo(:), dmosym(:)
    complex*16, allocatable :: itrfmog(:, :), itrfmou(:, :)
    logical                 :: write_itrfmo

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!! NOW MAKE FOCK MATRIX FOR IVO
!! fij = hij + SIGUMA_k (ij|kk)-(ik|kj)} i, j run over virtual spinors k runs occupied spinors except HOMO

    f = 0.0d+00

    if (rank == 0) print *, 'enter building fock matrix for IVO'

    if (nhomo == 0) then
        numh = count(ABS(caspt2_mo_energy(1:ninact + nact) - caspt2_mo_energy(nelec + ninact)) < 1.0d-01)
    else
        numh = nhomo
    end if

    if (rank == 0) print *, 'number of degeneracy of HOMO is', numh, DBLE(numh), 1.0d+00/DBLE(numh)

    do i = 1, nsec
        i0 = i + ninact + nact
        f(i, i) = caspt2_mo_energy(i0)
        do j = i, nsec
            j0 = j + ninact + nact
            do k = ninact + nact - numh + 1, ninact + nact

                if (k > ninact + nact - 2 .and. mod(nelec, 2) == 1) then
                    f(i, j) = f(i, j) - 0.5d+00*DCMPLX(int2r_f1(i0, j0, k, k), int2i_f1(i0, j0, k, k))/DBLE(numh)
                    f(i, j) = f(i, j) + 0.5d+00*DCMPLX(int2r_f2(i0, k, k, j0), int2i_f2(i0, k, k, j0))/DBLE(numh)
                else
                    f(i, j) = f(i, j) - DCMPLX(int2r_f1(i0, j0, k, k), int2i_f1(i0, j0, k, k))/DBLE(numh)
                    f(i, j) = f(i, j) + DCMPLX(int2r_f2(i0, k, k, j0), int2i_f2(i0, k, k, j0))/DBLE(numh)
                end if

            end do
        end do
    end do

    do i = 1, nsec
        do j = i, nsec
            f(j, i) = DCONJG(f(i, j))
        end do
    end do

    IMAX = nbas*lscom
    call open_formatted_file(unit=unit_dfpcmo, file='DFPCMO', status='old', optional_action='read')

! From DIRAC dirgp.F WRIPCMO (Write DHF-coefficients and eigenvalues )

    if (dirac_version == 21 .or. dirac_version == 22) then
        read (unit_dfpcmo, '(A150)') line0
    end if
    read (unit_dfpcmo, '(A150)') line1
    if (dirac_version == 21 .or. dirac_version == 22) then
        read (unit_dfpcmo, *) A, B, npg, neg, nbasg, npu, neu, nbasu
        if (rank == 0) print *, A, B, npg, neg, nbasg, npu, neu, nbasu
    else
        read (unit_dfpcmo, *) A, npg, neg, nbasg, npu, neu, nbasu
        if (rank == 0) print *, A, npg, neg, nbasg, npu, neu, nbasu
    end if
    read (unit_dfpcmo, '(A150)') line2

    if (rank == 0) print *, 'end reading information, symmetry information and energy'

    nsum = npg + neg + npu + neu

    ! ngu : the sum of gerade and ungerade
    ngu = (npg + neg)*nbasg + (npu + neu)*nbasu

    allocate (itrfmog(nbasg, neg - noccg - nvcutg))
    allocate (itrfmou(nbasu, neu - noccu - nvcutu))
    allocate (eval(nsum))
    allocate (syminfo(nsum))
    Allocate (BUF(ngu))

    itrfmog = 0.0d+00
    itrfmou = 0.0d+00
    BUF = 0.0d+00

    if (dirac_version == 21 .or. dirac_version == 22) then
        read (unit_dfpcmo, '(A150)') line3
    end if
    ! Read MO coefficient of DFPCMO
    write_itrfmo = .true.
    Do I = 1, ngu, 6
        Read (unit_dfpcmo, '(6F22.16)', iostat=iostat) BUF(I:I + 5)
        if (iostat /= 0) then
            write_itrfmo = .false.
            exit
        end if
    End do

    if (rank == 0) print *, 'end reading MO coefficient'
    if (write_itrfmo) then
        ! unoccupid, gerade, electron
        Do iao = 1, nbasg
            DO imo = 1, neg - noccg - nvcutg
                itrfmog(iao, imo) = BUF((npg + noccg + (imo - 1))*nbasg + iao)
            End do
        End do
        if (rank == 0) then
            call open_formatted_file(unit=unit_itrfmog, file='itrfmog_before', status='replace', optional_action='write')
            Do iao = 1, nbasg
                write (unit_itrfmog, *) (itrfmog(iao, imo), imo=1, neg - noccg - nvcutg)
            End do
            close (unit_itrfmog)
        end if
        ! unoccupid, ungerade, electron
        Do iao = 1, nbasu
            DO imo = 1, neu - noccu - nvcutu
                itrfmou(iao, imo) = BUF((npg + neg)*nbasg + (npu + noccu + (imo - 1))*nbasu + iao)
            End do
        End do
        if (rank == 0) then
            call open_formatted_file(unit=unit_itrfmou, file='itrfmou_before', status='replace', optional_action='write')
            Do iao = 1, nbasu
                write (unit_itrfmou, *) (itrfmou(iao, imo), imo=1, neu - noccu - nvcutu)
            End do
            close (unit_itrfmou)
        end if
    end if

    if (dirac_version == 21 .or. dirac_version == 22) then
        read (unit_dfpcmo, '(A150)') line4
    end if

    Read (unit_dfpcmo, *) eval
    if (rank == 0) then
        Do i = 1, nsum
            print *, "eval(i)", eval(i)
        End do
        print *, 'end reading eigenvalue'
    end if
    ! Read syminfo from DFPCMO
    if (dirac_version == 21 .or. dirac_version == 22) then
        read (unit_dfpcmo, '(A150)') line5
    end if

    Read (unit_dfpcmo, *) syminfo
    if (rank == 0) then
        Do i = 1, nsum
            print *, "syminfo(i)", syminfo(i)
        End do
        print *, 'end reading symmetry information2'
    end if
    close (unit_dfpcmo)

    if (rank == 0) then
        call open_formatted_file(unit=unit_buf, file='BUF_write', status='replace', optional_action='write')
        Do I = 1, ngu, 6
            Write (unit_buf, '(6F22.16)') BUF(I:I + 5)
        End do
        close (unit_buf)
    end if

! IVO calculation

! gerade

    Do isym = 1, nsymrpa, 2
        nv = count(irpmo(ninact + nact + 1:ninact + nact + nsec) == isym)

        Allocate (mosym(nv))
        Allocate (fsym(nv, nv))

        fsym = 0.0d+00
        nv = 0
        Do i = 1, nsec
            i0 = i + ninact + nact
            if (irpmo(i0) == isym) then
!                 if(irpamo(i0)==isym) then
                nv = nv + 1
                mosym(nv) = i
            end if
        end do

!C32h gerade
        if (isym <= nsymrpa/2) then
            nv0 = count(ABS(syminfo(npg + noccg + 1:npg + neg - nvcutg)) == isym)

            Allocate (dmosym(nv0))

            nv0 = 0
            Do i0 = npg + noccg + 1, npg + neg - nvcutg
                if (ABS(syminfo(i0)) == isym) then
!                 if(irpamo(i0)==isym) then
                    nv0 = nv0 + 1
                    dmosym(nv0) = i0
                end if
            end do

!C32h ungerade
        else
            nv0 = count(ABS(syminfo(npg + neg + npu + noccu + 1:npg + neg + npu + neu - nvcutu)) + nsymrpa/2 == isym)

            Allocate (dmosym(nv0))

            nv0 = 0
            Do i0 = npg + neg + npu + noccu + 1, npg + neg + npu + neu - nvcutu
                if (ABS(syminfo(i0)) + nsymrpa/2 == isym) then
!                 if(irpamo(i0)==isym) then
                    nv0 = nv0 + 1
                    dmosym(nv0) = i0
                end if
            end do
        end if

        Do i = 1, nv
            i0 = mosym(i)
            Do j = i, nv
                j0 = mosym(j)
                fsym(i, j) = f(i0, j0)
                fsym(j, i) = DCONJG(f(i0, j0))
!                    write(*,*)fsym(i,j)
            end do
        end do

        Allocate (wsym(nv))
        wsym = 0.0d+00
        cutoff = .FALSE.
        thresd = 0.0d+00

        call cdiag(fsym, nv, nv, wsym, thresd, cutoff)

        ! Gerade
        if (isym <= nsymrpa/2) then
            Allocate (coeff(nbasg, nv))
            Do i = 1, nv0
                i0 = dmosym(i) - npg - noccg
                coeff(:, i) = itrfmog(:, i0)
            End do

            ! Ungerade
        else
            Allocate (coeff(nbasu, nv))
            Do i = 1, nv0
                i0 = dmosym(i) - npg - neg - npu - noccu
                coeff(:, i) = itrfmou(:, i0)
            End do
        end if

        coeff(:, :) = MATMUL(coeff(:, :), fsym(:, :))

        ! Gerade
        if (isym <= nsymrpa/2) then
            Do i = 1, nv0
                i0 = dmosym(i) - npg - noccg
                itrfmog(:, i0) = coeff(:, i)
            End do

            ! Ungerade
        else
            Do i = 1, nv0
                i0 = dmosym(i) - npg - neg - npu - noccu
                itrfmou(:, i0) = coeff(:, i)
            End do
        end if

        deallocate (coeff)

        Do i = 1, nv
            i0 = mosym(i)
            if (rank == 0) print '(I4,F20.10)', i0, wsym(i)
        end do

        Do i = 1, nv
            i0 = mosym(i)
            if (rank == 0) then
                print *, ''
                print *, 'new ', i0 + ninact + nact, 'th ms consists of '
            end if
            Do j = 1, nv
                j0 = mosym(j)
                if (ABS(fsym(j, i))**2 > 1.0d-03) then
                    if (rank == 0) print '(I4,"  Weights ",F20.10)', j0 + ninact + nact, ABS(fsym(j, i))**2
                end if
            end do
        end do
        deallocate (fsym)
        deallocate (wsym)
        deallocate (mosym)
        deallocate (dmosym)
    end do

    ! unoccupid, gerade, electron
    Do iao = 1, nbasg
        DO imo = 1, neg - noccg - nvcutg
            BUF((npg + noccg + (imo - 1))*nbasg + iao) = DBLE(itrfmog(iao, imo))
        End do
    End do

    ! unoccupid, ungerade, electron
    Do iao = 1, nbasu
        DO imo = 1, neu - noccu - nvcutu
            BUF((npg + neg)*nbasg + (npu + noccu + (imo - 1))*nbasu + iao) = DBLE(itrfmou(iao, imo))
        End do
    End do

! Create new DFPCMO : DFPCMONEW
    if (rank == 0) then
        call open_formatted_file(unit=unit_dfpcmo, file='DFPCMONEW', status='replace', optional_action="write")
        if (dirac_version == 21 .or. dirac_version == 22) then
            write (unit_dfpcmo, '(A150)') line0
        end if
        write (unit_dfpcmo, '(A150)') line1
        if (dirac_version == 21 .or. dirac_version == 22) then
            write (unit_dfpcmo, '(8(X,I0))') A, B, npg, neg, nbasg, npu, neu, nbasu
        else
            write (unit_dfpcmo, '(7(X,I0))') A, npg, neg, nbasg, npu, neu, nbasu
        end if
        write (unit_dfpcmo, '(A150)') line2
        if (dirac_version == 21 .or. dirac_version == 22) then
            write (unit_dfpcmo, '(A150)') line3
        end if
        Do I = 1, ngu, 6
            Write (unit_dfpcmo, '(6F22.16)') BUF(I:I + 5)
        End do
        if (dirac_version == 21 .or. dirac_version == 22) then
            write (unit_dfpcmo, '(A150)') line4
        end if

        write (unit_dfpcmo, '(6E22.12)') eval
        if (dirac_version == 21 .or. dirac_version == 22) then
            write (unit_dfpcmo, '(A150)') line5
        end if
        write (unit_dfpcmo, '(66(X,I0))') (syminfo(i), i=1, nsum)
        close (unit_dfpcmo)
    end if
    if (rank == 0) print *, 'fockivo_co end'
    deallocate (itrfmog, itrfmou)
    deallocate (BUF)
    deallocate (eval, syminfo)
end subroutine fockivo_co
