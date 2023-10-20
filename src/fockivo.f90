module fockivo

    implicit none

    private
    public fockivo_cmplx, fockivo_real

contains

    subroutine fockivo_real
        use module_global_variables
        use module_file_manager
        use module_index_utils, only: convert_secondary_to_global_idx

        Implicit none

        integer      :: j, i, k, i0, j0
        integer      :: isym, nv, numh
        integer      :: imo, iao, iostat
        integer      :: unit_dfpcmo, unit_buf
        real(8)       :: thresd
        real(8), allocatable  :: fsym(:, :) ! Symmetrized fock_ivo_matrix for particular irrep
        real(8), allocatable  :: coeff(:, :), itrfmo(:, :)
        real(8), allocatable      :: BUF(:)  ! One dimensional array representing MO coeff. read from DFPCMO
        real(8), allocatable      :: wsym(:), eval(:)
        integer, allocatable     :: mosym(:)
        character*150 :: line0, line1, line2, line3, line4, line5, format_str

        ! for new code of IVO
        integer :: total_ao, total_mo
        integer :: nv0, A, B ! A and B are dammy indices written in DFPCMO, A is nfsym in DIRAC
        integer :: idx_irrep, start_isym, end_isym
        integer, allocatable :: syminfo(:), dmosym(:)
        ! complex*16, allocatable :: itrfmo(:, :)
        logical                 :: write_itrfmo, is_all_syminfo_zero
        integer :: juck_up_idx, num_ao, num_mo, num_virtual_mo
        integer :: mo_start_idx, mo_end_idx, isym_for_syminfo
        integer :: positronic_mo(2), electronic_mo(2), basis_ao(2), basis_all(2), mo(2)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!! NOW MAKE FOCK MATRIX FOR IVO
!! fij = hij + SIGUMA_k (ij|kk)-(ik|kj)} i, j run over virtual spinors k runs occupied spinors except HOMO

        if (allocated(fock_real)) then
            call memminus(kind(fock_real), size(fock_real), 1)
            deallocate (fock_real)
        end if
        allocate (fock_real(nsec, nsec)); call memplus(kind(fock_real), size(fock_real), 1)

        fock_real(:, :) = 0.0d+00; positronic_mo(:) = 0; electronic_mo(:) = 0; basis_ao(:) = 0; basis_all(:) = 0; mo(:) = 0

        if (rank == 0) print *, 'enter building fock matrix for IVO'

        if (nhomo == 0) then
            numh = count(ABS(caspt2_mo_energy(1:global_act_end) - caspt2_mo_energy(nelec + ninact)) < 1.0d-01)
        else
            numh = nhomo
        end if

        if (rank == 0) print *, 'number of degeneracy of HOMO is', numh, DBLE(numh), 1.0d+00/DBLE(numh)

        ! Create Fock matrix (only virtual)
        do i = 1, nsec
            i0 = convert_secondary_to_global_idx(i)
            fock_real(i, i) = caspt2_mo_energy(i0)
            do j = i, nsec
                j0 = convert_secondary_to_global_idx(j)
                do k = global_act_end - numh + 1, global_act_end

                    if (k > global_act_end - 2 .and. mod(nelec, 2) == 1) then
                        fock_real(i, j) = fock_real(i, j) - 0.5d+00*int2r_f1(i0, j0, k, k)/DBLE(numh)
                        fock_real(i, j) = fock_real(i, j) + 0.5d+00*int2r_f2(i0, k, k, j0)/DBLE(numh)
                    else
                        fock_real(i, j) = fock_real(i, j) + int2r_f2(i0, k, k, j0)/DBLE(numh)
                        fock_real(i, j) = fock_real(i, j) - int2r_f1(i0, j0, k, k)/DBLE(numh)
                    end if

                end do
            end do
        end do

        call open_formatted_file(unit=unit_dfpcmo, file='DFPCMO', status='old', optional_action='read')
        rewind (unit_dfpcmo)

! From DIRAC dirgp.F WRIPCMO (Write DHF-coefficients and eigenvalues )

        if (dirac_version >= 21) then
            read (unit_dfpcmo, '(A150)') line0
        end if
        read (unit_dfpcmo, '(A150)') line1
        if (dirac_version >= 21) then
            ! A is nfsym2 in DIRAC (https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirgp.F#L77-78)
            ! A is 1 or 2
            ! If B is 1, COEFS in DFPCMO is written as Complex number, if B is 2, COEFS in DFPCMO is written as Real number
            read (unit_dfpcmo, *) A, B, (positronic_mo(idx_irrep), electronic_mo(idx_irrep), basis_ao(idx_irrep), idx_irrep=1, A)
        else
            read (unit_dfpcmo, *) A, (positronic_mo(idx_irrep), electronic_mo(idx_irrep), basis_ao(idx_irrep), idx_irrep=1, A)
        end if
        read (unit_dfpcmo, '(A150)') line2

        if (rank == 0) print *, 'end reading information, symmetry information and energy'

        basis_all = (positronic_mo + electronic_mo)*basis_ao
        mo = positronic_mo + electronic_mo
        total_mo = sum(positronic_mo) + sum(electronic_mo)
        total_ao = sum(basis_all)
        allocate (eval(total_mo))
        allocate (syminfo(total_mo))
        Allocate (BUF(total_ao))

        BUF = 0.0d+00
        if (dirac_version >= 21) then
            read (unit_dfpcmo, '(A150)') line3
        end if
        ! Read MO coefficient of DFPCMO
        write_itrfmo = .true.
        read (unit_dfpcmo, *, iostat=iostat) BUF

        if (rank == 0) print *, 'end reading MO coefficient'

        if (dirac_version >= 21) then
            read (unit_dfpcmo, '(A150)') line4
        end if

        Read (unit_dfpcmo, *) eval
        if (rank == 0) then
            Do i = 1, total_mo
                print *, "eval(i)", eval(i)
            End do
            print *, 'end reading eigenvalue'
        end if
        if (dirac_version >= 21) then
            read (unit_dfpcmo, '(A150)') line5
        end if

        Read (unit_dfpcmo, *) syminfo
        if (rank == 0) then
            Do i = 1, total_mo
                print *, "syminfo(i)", syminfo(i)
            End do
            print *, 'end reading symmetry information2'
        end if
        close (unit_dfpcmo)

        if (rank == 0) then
            call open_formatted_file(unit=unit_buf, file='BUF_write', status='replace', optional_action='write')
            Do I = 1, total_ao, 6
                Write (unit_buf, '(6F22.16)') BUF(I:I + 5)
            End do
            close (unit_buf)
        end if

! IVO calculation (C1 symmetry is not supported)
        do idx_irrep = 1, A
            num_ao = basis_ao(idx_irrep)
            if (idx_irrep == 1) then
                num_virtual_mo = electronic_mo(idx_irrep) - occ_mo_num(idx_irrep) - vcut_mo_num(idx_irrep)
                start_isym = 1
                end_isym = nsymrpa/2
                juck_up_idx = (positronic_mo(idx_irrep) + occ_mo_num(idx_irrep))*basis_ao(idx_irrep)
                num_mo = electronic_mo(idx_irrep) - occ_mo_num(idx_irrep) - vcut_mo_num(idx_irrep)
            else
                num_virtual_mo = electronic_mo(idx_irrep) - occ_mo_num(idx_irrep) - vcut_mo_num(idx_irrep)
                start_isym = nsymrpa/2 + 1
                end_isym = nsymrpa
                juck_up_idx = basis_all(1) + (positronic_mo(idx_irrep) + occ_mo_num(idx_irrep))*basis_ao(idx_irrep)
                num_mo = electronic_mo(idx_irrep) - occ_mo_num(idx_irrep) - vcut_mo_num(idx_irrep)
            end if

            allocate (itrfmo(basis_ao(idx_irrep), num_virtual_mo))
            itrfmo(:, :) = 0.0d+00
            call create_itrfmo

            num_mo = num_virtual_mo

            Do isym = start_isym, end_isym, 2
                nv = count(irpamo(global_sec_start:global_sec_end) == isym)

                Allocate (mosym(nv))
                Allocate (fsym(nv, nv))

                fsym = 0.0d+00
                nv = 0
                Do i = 1, nsec
                    i0 = convert_secondary_to_global_idx(i)
                    if (irpamo(i0) == isym) then
                        nv = nv + 1
                        mosym(nv) = i
                    end if
                end do

                if (idx_irrep == 1) then
                    mo_start_idx = positronic_mo(idx_irrep) + occ_mo_num(idx_irrep) + 1
                    mo_end_idx = positronic_mo(idx_irrep) + electronic_mo(idx_irrep) - vcut_mo_num(idx_irrep)
                    isym_for_syminfo = isym
                else
                    mo_start_idx = mo(1) + positronic_mo(idx_irrep) + occ_mo_num(idx_irrep) + 1
                    mo_end_idx = mo(1) + positronic_mo(idx_irrep) + electronic_mo(idx_irrep) - vcut_mo_num(idx_irrep)
                    isym_for_syminfo = isym - nsymrpa/2
                end if

                if (all(syminfo(mo_start_idx:mo_end_idx) == 0)) then
                    print *, 'all syminfo is zero, idx_irrep = ', idx_irrep
                    nv0 = mo_end_idx - mo_start_idx + 1
                    is_all_syminfo_zero = .true.
                else
                    nv0 = count(ABS(syminfo(mo_start_idx:mo_end_idx)) == isym_for_syminfo)
                    is_all_syminfo_zero = .false.
                end if
                Allocate (dmosym(nv0))
                call create_dmosym(mo_start_idx, mo_end_idx)

                Do i = 1, nv
                    i0 = mosym(i)
                    Do j = i, nv
                        j0 = mosym(j)
                        fsym(i, j) = fock_real(i0, j0)
                    end do
                end do

                Allocate (wsym(nv))
                wsym = 0.0d+00
                thresd = 0.0d+00

                call rdiag(fsym, nv, nv, wsym, thresd)

                ! Gerade
                allocate (coeff(basis_ao(idx_irrep), nv))
                coeff(:, :) = 0.0d+00
                if (idx_irrep == 1) then
                    juck_up_idx = positronic_mo(idx_irrep) + occ_mo_num(idx_irrep)
                else
                    juck_up_idx = mo(1) + positronic_mo(idx_irrep) + occ_mo_num(idx_irrep)
                end if
                call get_coeff

                coeff(:, :) = MATMUL(coeff(:, :), fsym(:, :))
                if (idx_irrep == 1) then
                    juck_up_idx = positronic_mo(idx_irrep) + occ_mo_num(idx_irrep)
                else
                    juck_up_idx = mo(1) + positronic_mo(idx_irrep) + occ_mo_num(idx_irrep)
                end if

                call write_back_itrfmo

                deallocate (coeff)

                Do i = 1, nv
                    i0 = mosym(i)
                    if (rank == 0) print '(I4,F20.10)', i0, wsym(i)
                end do

                Do i = 1, nv
                    i0 = mosym(i)
                    if (rank == 0) then
                        print *, ''
                        print *, 'new ', convert_secondary_to_global_idx(i0), 'th ms consists of '
                    end if
                    Do j = 1, nv
                        j0 = mosym(j)
                        if (ABS(fsym(j, i))**2 > 1.0d-03) then
                            if (rank == 0) print '(I4,"  Weights ",F20.10)', convert_secondary_to_global_idx(j0), ABS(fsym(j, i))**2
                        end if
                    end do
                end do
                deallocate (fsym)
                deallocate (wsym)
                deallocate (mosym)
                deallocate (dmosym)
            end do
            num_ao = basis_ao(idx_irrep)
            if (idx_irrep == 1) then
                juck_up_idx = (positronic_mo(idx_irrep) + occ_mo_num(idx_irrep))*basis_ao(idx_irrep)
                num_mo = electronic_mo(idx_irrep) - occ_mo_num(idx_irrep) - vcut_mo_num(idx_irrep)
            else
                juck_up_idx = basis_all(1) + (positronic_mo(idx_irrep) + occ_mo_num(idx_irrep))*basis_ao(idx_irrep)
                num_mo = electronic_mo(idx_irrep) - occ_mo_num(idx_irrep) - vcut_mo_num(idx_irrep)
            end if

            do iao = 1, num_ao
                do imo = 1, num_mo
                    BUF(juck_up_idx + (imo - 1)*num_ao + iao) = DBLE(itrfmo(iao, imo))
                end do
            end do
            deallocate (itrfmo)
        end do

! Create new DFPCMO : DFPCMONEW
        if (rank == 0) then
            call open_formatted_file(unit=unit_dfpcmo, file='DFPCMONEW', status='replace', optional_action="write")
            if (dirac_version >= 21) then
                write (unit_dfpcmo, '(A150)') line0
            end if
            write (unit_dfpcmo, '(A150)') line1
            if (dirac_version >= 21) then
                if (A == 1) then
                    format_str = '(5(X,I0))'
                else
                    format_str = '(8(X,I0))'
                end if
                write (unit_dfpcmo, format_str) A, B, (positronic_mo(i), electronic_mo(i), basis_ao(i), i=1, A)
            else
                if (A == 1) then
                    format_str = '(4(X,I0))'
                else
                    format_str = '(7(X,I0))'
                end if
                write (unit_dfpcmo, format_str) A, (positronic_mo(i), electronic_mo(i), basis_ao(i), i=1, A)
            end if
            write (unit_dfpcmo, '(A150)') line2
            if (dirac_version >= 21) then
                write (unit_dfpcmo, '(A150)') line3
            end if
            Do I = 1, total_ao, 6
                Write (unit_dfpcmo, '(6F22.16)') BUF(I:I + 5)
            End do
            if (dirac_version >= 21) then
                write (unit_dfpcmo, '(A150)') line4
            end if

            write (unit_dfpcmo, '(6E22.12)') eval
            if (dirac_version >= 21) then
                write (unit_dfpcmo, '(A150)') line5
            end if
            write (unit_dfpcmo, '(66(X,I0))') (syminfo(i), i=1, total_mo)
            close (unit_dfpcmo)
        end if
        if (rank == 0) print *, 'fockivo end'
        deallocate (BUF)
        deallocate (eval, syminfo)
        call memminus(kind(fock_real), size(fock_real), 1); deallocate (fock_real)
    contains
        subroutine create_itrfmo
            implicit none

            do iao = 1, num_ao
                do imo = 1, num_mo
                    itrfmo(iao, imo) = BUF(juck_up_idx + (imo - 1)*num_ao + iao)
                end do
            end do
        end subroutine create_itrfmo

        subroutine create_dmosym(start_idx, end_idx)
            use module_error, only: stop_with_errorcode
            implicit none
            integer, intent(in) :: start_idx, end_idx
            integer :: idx, cnt

            dmosym(:) = 0
            cnt = 0
            do idx = start_idx, end_idx
                if (is_all_syminfo_zero .or. abs(syminfo(idx)) == isym_for_syminfo) then
                    cnt = cnt + 1
                    dmosym(cnt) = idx
                end if
            end do

            ! Validate
            if (cnt /= size(dmosym)) then
                if (rank == 0) print '(a,i0,a,i0)', &
                    'Error in create_dmosym, cnt /= size(dmosym). cnt = ', cnt, 'size(dmosym) = ', size(dmosym)
                call stop_with_errorcode(1)
            end if
        end subroutine create_dmosym

        subroutine get_coeff
            implicit none

            Do i = 1, nv0
                i0 = dmosym(i) - juck_up_idx
                coeff(:, i) = itrfmo(:, i0)
            End do
        end subroutine get_coeff

        subroutine write_back_itrfmo
            implicit none

            Do i = 1, nv0
                i0 = dmosym(i) - juck_up_idx
                itrfmo(:, i0) = coeff(:, i)
            End do
        end subroutine write_back_itrfmo
    end subroutine fockivo_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE fockivo_cmplx ! TO MAKE FOCK MATRIX for IVO

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_file_manager
        use module_index_utils, only: convert_secondary_to_global_idx

        Implicit NONE

        integer      :: j, i, k, i0, j0
        integer      :: isym, nv, numh
        integer      :: imo, iao, iostat
        integer      :: unit_dfpcmo, unit_buf
        real(8)       :: thresd
        complex*16, allocatable  :: fsym(:, :) ! Symmetrized fock_ivo_matrix for particular irrep
        complex*16, allocatable  :: coeff(:, :), itrfmo(:, :)
        complex*16, allocatable      :: BUF(:)  ! One dimensional array representing MO coeff. read from DFPCMO
        real(8), allocatable      :: wsym(:), eval(:), BUF_READWRITE(:)
        integer, allocatable     :: mosym(:)
        character*150 :: line0, line1, line2, line3, line4, line5, format_str

        ! for new code of IVO
        integer :: total_ao, total_mo
        integer :: nv0, A, B ! A and B are dammy indices written in DFPCMO, A is nfsym in DIRAC
        integer :: idx_irrep, start_isym, end_isym
        integer, allocatable :: syminfo(:), dmosym(:)
        ! complex*16, allocatable :: itrfmo(:, :)
        logical                 :: write_itrfmo, is_all_syminfo_zero
        integer :: juck_up_idx, num_ao, num_mo, num_virtual_mo
        integer :: mo_start_idx, mo_end_idx, isym_for_syminfo
        integer :: positronic_mo(2), electronic_mo(2), basis_ao(2), basis_all(2), mo(2)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!! NOW MAKE FOCK MATRIX FOR IVO
!! fij = hij + SIGUMA_k (ij|kk)-(ik|kj)} i, j run over virtual spinors k runs occupied spinors except HOMO
        if (allocated(fock_cmplx)) then
            call memminus(KIND(fock_cmplx), SIZE(fock_cmplx), 2)
            deallocate (fock_cmplx)
        end if
        Allocate (fock_cmplx(nsec, nsec)); Call memplus(KIND(fock_cmplx), SIZE(fock_cmplx), 2)
        fock_cmplx(:, :) = 0.0d+00; positronic_mo(:) = 0; electronic_mo(:) = 0; basis_ao(:) = 0; basis_all(:) = 0; mo(:) = 0

        if (rank == 0) print *, 'enter building fock matrix for IVO'

        if (nhomo == 0) then
            numh = count(ABS(caspt2_mo_energy(1:global_act_end) - caspt2_mo_energy(nelec + ninact)) < 1.0d-01)
        else
            numh = nhomo
        end if

        if (rank == 0) print *, 'number of degeneracy of HOMO is', numh, DBLE(numh), 1.0d+00/DBLE(numh)

        ! Create Fock matrix (only virtual)
        do i = 1, nsec
            i0 = convert_secondary_to_global_idx(i)
            fock_cmplx(i, i) = caspt2_mo_energy(i0)
            do j = i, nsec
                j0 = convert_secondary_to_global_idx(j)
                do k = global_act_end - numh + 1, global_act_end

                    if (k > global_act_end - 2 .and. mod(nelec, 2) == 1) then
                        fock_cmplx(i, j) = fock_cmplx(i, j) - &
                                           0.5d+00*DCMPLX(int2r_f1(i0, j0, k, k), int2i_f1(i0, j0, k, k))/DBLE(numh)
                        fock_cmplx(i, j) = fock_cmplx(i, j) + &
                                           0.5d+00*DCMPLX(int2r_f2(i0, k, k, j0), int2i_f2(i0, k, k, j0))/DBLE(numh)
                    else
                        fock_cmplx(i, j) = fock_cmplx(i, j) - DCMPLX(int2r_f1(i0, j0, k, k), int2i_f1(i0, j0, k, k))/DBLE(numh)
                        fock_cmplx(i, j) = fock_cmplx(i, j) + DCMPLX(int2r_f2(i0, k, k, j0), int2i_f2(i0, k, k, j0))/DBLE(numh)
                    end if

                end do
            end do
        end do
        ! Take conjugate
        do i = 1, nsec
            do j = i, nsec
                fock_cmplx(j, i) = DCONJG(fock_cmplx(i, j))
            end do
        end do

        call open_formatted_file(unit=unit_dfpcmo, file='DFPCMO', status='old', optional_action='read')
        rewind (unit_dfpcmo)

! From DIRAC dirgp.F WRIPCMO (Write DHF-coefficients and eigenvalues )

        if (dirac_version >= 21) then
            read (unit_dfpcmo, '(A150)') line0
        end if
        read (unit_dfpcmo, '(A150)') line1
        if (dirac_version >= 21) then
            ! A is nfsym2 in DIRAC (https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirgp.F#L77-78)
            ! A is 1 or 2
            ! If B is 1, COEFS in DFPCMO is written as Complex number, if B is 2, COEFS in DFPCMO is written as Real number
            read (unit_dfpcmo, *) A, B, (positronic_mo(idx_irrep), electronic_mo(idx_irrep), basis_ao(idx_irrep), idx_irrep=1, A)
        else
            read (unit_dfpcmo, *) A, (positronic_mo(idx_irrep), electronic_mo(idx_irrep), basis_ao(idx_irrep), idx_irrep=1, A)
        end if
        read (unit_dfpcmo, '(A150)') line2

        if (rank == 0) print *, 'end reading information, symmetry information and energy'

        basis_all = (positronic_mo + electronic_mo)*basis_ao
        mo = positronic_mo + electronic_mo
        total_mo = sum(positronic_mo) + sum(electronic_mo)
        total_ao = sum(basis_all)
        allocate (eval(total_mo))
        allocate (syminfo(total_mo))
        Allocate (BUF(total_ao))

        BUF = 0.0d+00
        if (dirac_version >= 21) then
            read (unit_dfpcmo, '(A150)') line3
        end if
        ! Read MO coefficient of DFPCMO
        write_itrfmo = .true.
        print *, 'total_ao = ', total_ao
        allocate (BUF_READWRITE(total_ao*2))
        read (unit_dfpcmo, *, iostat=iostat) BUF_READWRITE
        do i = 1, total_ao
            BUF(i) = DCMPLX(BUF_READWRITE(i), BUF_READWRITE(total_ao + i))
        end do

        if (rank == 0) then
            print *, 'end reading MO coefficient'
            call open_formatted_file(unit=unit_buf, file='BUF_write', status='replace', optional_action='write')
            Do I = 1, total_ao, 6
                Write (unit_buf, '(6F22.16)') BUF_READWRITE(I:I + 5)
            End do
            close (unit_buf)
        end if

        if (dirac_version >= 21) then
            read (unit_dfpcmo, '(A150)') line4
        end if

        Read (unit_dfpcmo, *) eval
        if (rank == 0) then
            Do i = 1, total_mo
                print *, "eval(i)", eval(i)
            End do
            print *, 'end reading eigenvalue'
        end if
        if (dirac_version >= 21) then
            read (unit_dfpcmo, '(A150)') line5
        end if

        Read (unit_dfpcmo, *) syminfo
        if (rank == 0) then
            Do i = 1, total_mo
                print *, "syminfo(i)", syminfo(i)
            End do
            print *, 'end reading symmetry information2'
        end if
        close (unit_dfpcmo)

! IVO calculation (C1 symmetry is not supported)
        do idx_irrep = 1, A
            num_ao = basis_ao(idx_irrep)
            if (idx_irrep == 1) then
                num_virtual_mo = electronic_mo(idx_irrep) - occ_mo_num(idx_irrep) - vcut_mo_num(idx_irrep)
                start_isym = 1
                end_isym = nsymrpa/2
                juck_up_idx = (positronic_mo(idx_irrep) + occ_mo_num(idx_irrep))*basis_ao(idx_irrep)
                num_mo = electronic_mo(idx_irrep) - occ_mo_num(idx_irrep) - vcut_mo_num(idx_irrep)
            else
                num_virtual_mo = electronic_mo(idx_irrep) - occ_mo_num(idx_irrep) - vcut_mo_num(idx_irrep)
                start_isym = nsymrpa/2 + 1
                end_isym = nsymrpa
                juck_up_idx = basis_all(1) + (positronic_mo(idx_irrep) + occ_mo_num(idx_irrep))*basis_ao(idx_irrep)
                num_mo = electronic_mo(idx_irrep) - occ_mo_num(idx_irrep) - vcut_mo_num(idx_irrep)
            end if

            allocate (itrfmo(basis_ao(idx_irrep), num_virtual_mo))
            itrfmo(:, :) = 0.0d+00
            call create_itrfmo

            num_mo = num_virtual_mo

            Do isym = start_isym, end_isym, 2
                nv = count(irpamo(global_sec_start:global_sec_end) == isym)

                Allocate (mosym(nv))
                Allocate (fsym(nv, nv))

                fsym = 0.0d+00
                nv = 0
                Do i = 1, nsec
                    i0 = convert_secondary_to_global_idx(i)
                    if (irpamo(i0) == isym) then
                        nv = nv + 1
                        mosym(nv) = i
                    end if
                end do

                if (idx_irrep == 1) then
                    mo_start_idx = positronic_mo(idx_irrep) + occ_mo_num(idx_irrep) + 1
                    mo_end_idx = positronic_mo(idx_irrep) + electronic_mo(idx_irrep) - vcut_mo_num(idx_irrep)
                    isym_for_syminfo = isym
                else
                    mo_start_idx = mo(1) + positronic_mo(idx_irrep) + occ_mo_num(idx_irrep) + 1
                    mo_end_idx = mo(1) + positronic_mo(idx_irrep) + electronic_mo(idx_irrep) - vcut_mo_num(idx_irrep)
                    isym_for_syminfo = isym - nsymrpa/2
                end if

                if (all(syminfo(mo_start_idx:mo_end_idx) == 0)) then
                    print *, 'all syminfo is zero, idx_irrep = ', idx_irrep
                    nv0 = mo_end_idx - mo_start_idx + 1
                    is_all_syminfo_zero = .true.
                else
                    nv0 = count(ABS(syminfo(mo_start_idx:mo_end_idx)) == isym_for_syminfo)
                    is_all_syminfo_zero = .false.
                end if
                Allocate (dmosym(nv0))
                call create_dmosym(mo_start_idx, mo_end_idx)

                Do i = 1, nv
                    i0 = mosym(i)
                    Do j = i, nv
                        j0 = mosym(j)
                        fsym(i, j) = fock_cmplx(i0, j0)
                        fsym(j, i) = DCONJG(fock_cmplx(i0, j0))
                    end do
                end do

                Allocate (wsym(nv))
                wsym = 0.0d+00
                thresd = 0.0d+00

                call cdiag(fsym, nv, nv, wsym, thresd)

                ! Gerade
                allocate (coeff(basis_ao(idx_irrep), nv))
                coeff(:, :) = 0.0d+00
                if (idx_irrep == 1) then
                    juck_up_idx = positronic_mo(idx_irrep) + occ_mo_num(idx_irrep)
                else
                    juck_up_idx = mo(1) + positronic_mo(idx_irrep) + occ_mo_num(idx_irrep)
                end if
                call get_coeff

                coeff(:, :) = MATMUL(coeff(:, :), fsym(:, :))
                if (idx_irrep == 1) then
                    juck_up_idx = positronic_mo(idx_irrep) + occ_mo_num(idx_irrep)
                else
                    juck_up_idx = mo(1) + positronic_mo(idx_irrep) + occ_mo_num(idx_irrep)
                end if

                call write_back_itrfmo

                deallocate (coeff)

                Do i = 1, nv
                    i0 = mosym(i)
                    if (rank == 0) print '(I4,F20.10)', i0, wsym(i)
                end do

                Do i = 1, nv
                    i0 = mosym(i)
                    if (rank == 0) then
                        print *, ''
                        print *, 'new ', convert_secondary_to_global_idx(i0), 'th ms consists of '
                    end if
                    Do j = 1, nv
                        j0 = mosym(j)
                        if (ABS(fsym(j, i))**2 > 1.0d-03) then
                            if (rank == 0) print '(I4,"  Weights ",F20.10)', convert_secondary_to_global_idx(j0), ABS(fsym(j, i))**2
                        end if
                    end do
                end do
                deallocate (fsym)
                deallocate (wsym)
                deallocate (mosym)
                deallocate (dmosym)
            end do
            num_ao = basis_ao(idx_irrep)
            if (idx_irrep == 1) then
                juck_up_idx = (positronic_mo(idx_irrep) + occ_mo_num(idx_irrep))*basis_ao(idx_irrep)
                num_mo = electronic_mo(idx_irrep) - occ_mo_num(idx_irrep) - vcut_mo_num(idx_irrep)
            else
                juck_up_idx = basis_all(1) + (positronic_mo(idx_irrep) + occ_mo_num(idx_irrep))*basis_ao(idx_irrep)
                num_mo = electronic_mo(idx_irrep) - occ_mo_num(idx_irrep) - vcut_mo_num(idx_irrep)
            end if

            do iao = 1, num_ao
                do imo = 1, num_mo
                    BUF(juck_up_idx + (imo - 1)*num_ao + iao) = DBLE(itrfmo(iao, imo))
                end do
            end do
            deallocate (itrfmo)
        end do

! Create new DFPCMO : DFPCMONEW
        if (rank == 0) then
            call open_formatted_file(unit=unit_dfpcmo, file='DFPCMONEW', status='replace', optional_action="write")
            if (dirac_version >= 21) then
                write (unit_dfpcmo, '(A150)') line0
            end if
            write (unit_dfpcmo, '(A150)') line1
            if (dirac_version >= 21) then
                if (A == 1) then
                    format_str = '(5(X,I0))'
                else
                    format_str = '(8(X,I0))'
                end if
                write (unit_dfpcmo, format_str) A, B, (positronic_mo(i), electronic_mo(i), basis_ao(i), i=1, A)
            else
                if (A == 1) then
                    format_str = '(4(X,I0))'
                else
                    format_str = '(7(X,I0))'
                end if
                write (unit_dfpcmo, format_str) A, (positronic_mo(i), electronic_mo(i), basis_ao(i), i=1, A)
            end if
            write (unit_dfpcmo, '(A150)') line2
            if (dirac_version >= 21) then
                write (unit_dfpcmo, '(A150)') line3
            end if
            do i = 1, total_ao
                BUF_READWRITE(i) = real(BUF(i), kind=kind(BUF))
                BUF_READWRITE(total_ao + i) = aimag(BUF(i))
            end do
            Do I = 1, total_ao, 6
                Write (unit_dfpcmo, '(6F22.16)') BUF_READWRITE(I:I + 5)
            End do
            if (dirac_version >= 21) then
                write (unit_dfpcmo, '(A150)') line4
            end if

            write (unit_dfpcmo, '(6E22.12)') eval
            if (dirac_version >= 21) then
                write (unit_dfpcmo, '(A150)') line5
            end if
            write (unit_dfpcmo, '(66(X,I0))') (syminfo(i), i=1, total_mo)
            close (unit_dfpcmo)
        end if
        if (rank == 0) print *, 'fockivo end'
        deallocate (BUF)
        deallocate (eval, syminfo)
        deallocate (BUF_READWRITE)
        Call memminus(KIND(fock_cmplx), SIZE(fock_cmplx), 2); deallocate (fock_cmplx)
    contains
        subroutine create_itrfmo
            implicit none

            do iao = 1, num_ao
                do imo = 1, num_mo
                    itrfmo(iao, imo) = BUF(juck_up_idx + (imo - 1)*num_ao + iao)
                end do
            end do
        end subroutine create_itrfmo

        subroutine create_dmosym(start_idx, end_idx)
            use module_error, only: stop_with_errorcode
            implicit none
            integer, intent(in) :: start_idx, end_idx
            integer :: idx, cnt

            dmosym(:) = 0
            cnt = 0
            do idx = start_idx, end_idx
                if (is_all_syminfo_zero .or. abs(syminfo(idx)) == isym_for_syminfo) then
                    cnt = cnt + 1
                    dmosym(cnt) = idx
                end if
            end do

            ! Validate
            if (cnt /= size(dmosym)) then
                if (rank == 0) print '(a,i0,a,i0)', &
                    'Error in create_dmosym, cnt /= size(dmosym). cnt = ', cnt, 'size(dmosym) = ', size(dmosym)
                call stop_with_errorcode(1)
            end if
        end subroutine create_dmosym

        subroutine get_coeff
            implicit none

            Do i = 1, nv0
                i0 = dmosym(i) - juck_up_idx
                coeff(:, i) = itrfmo(:, i0)
            End do
        end subroutine get_coeff

        subroutine write_back_itrfmo
            implicit none

            Do i = 1, nv0
                i0 = dmosym(i) - juck_up_idx
                itrfmo(:, i0) = coeff(:, i)
            End do
        end subroutine write_back_itrfmo
    end subroutine fockivo_cmplx
end module fockivo
