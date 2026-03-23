! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE fockivo ! TO MAKE FOCK MATRIX for IVO

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_cmo_variables
    use module_cmo_handler, only: ivo_cmo_read, ivo_cmo_write
    use module_global_variables
    use module_file_manager
    use module_index_utils, only: convert_secondary_to_global_idx
    use module_realonly, only: realonly

    Implicit NONE

    integer      :: j, i, k, i0, j0
    integer      :: isym, nv, numh
    integer      :: imo, iao, unit_buf
    real(8)      :: thresd
    complex*16, allocatable  :: fsym(:, :) ! Symmetrized fock_ivo_matrix for particular irrep
    complex*16, allocatable  :: coeff(:, :)
    real(8), allocatable     :: wsym(:)
    integer, allocatable     :: mosym(:)

    integer :: nv0, idx_irrep, start_isym, end_isym
    logical :: is_all_syminfo_zero
    integer :: juck_up_idx, num_ao, num_mo, num_virtual_mo
    integer :: mo_start_idx, mo_end_idx, isym_for_syminfo
    integer :: kp, mj, ll, indi, idx
    integer, allocatable :: dmosym(:)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!! NOW MAKE FOCK MATRIX FOR IVO
!! fij = hij + SIGUMA_k (ij|kk)-(ik|kj)} i, j run over virtual spinors k runs occupied spinors except HOMO

    fock_cmplx = 0.0d+00
    positronic_mo(:) = 0; electronic_mo(:) = 0; basis_ao(:) = 0; basis_all(:) = 0; mo(:) = 0

    if (debug .and. rank == 0) print *, 'enter building fock matrix for IVO'

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
                    if (realonly%is_realonly()) then
                        fock_cmplx(i, j) = fock_cmplx(i, j) - 0.5d+00*int2r_f1(i0, j0, k, k)/DBLE(numh)
                        fock_cmplx(i, j) = fock_cmplx(i, j) + 0.5d+00*int2r_f2(i0, k, k, j0)/DBLE(numh)
                    else
                        fock_cmplx(i, j) = fock_cmplx(i, j) - &
                                           0.5d+00*DCMPLX(int2r_f1(i0, j0, k, k), int2i_f1(i0, j0, k, k))/DBLE(numh)
                        fock_cmplx(i, j) = fock_cmplx(i, j) + &
                                           0.5d+00*DCMPLX(int2r_f2(i0, k, k, j0), int2i_f2(i0, k, k, j0))/DBLE(numh)
                    end if
                else
                    if (realonly%is_realonly()) then
                        fock_cmplx(i, j) = fock_cmplx(i, j) + int2r_f2(i0, k, k, j0)/DBLE(numh)
                        fock_cmplx(i, j) = fock_cmplx(i, j) - int2r_f1(i0, j0, k, k)/DBLE(numh)
                    else
                        fock_cmplx(i, j) = fock_cmplx(i, j) - DCMPLX(int2r_f1(i0, j0, k, k), int2i_f1(i0, j0, k, k))/DBLE(numh)
                        fock_cmplx(i, j) = fock_cmplx(i, j) + DCMPLX(int2r_f2(i0, k, k, j0), int2i_f2(i0, k, k, j0))/DBLE(numh)
                    end if
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

    call reset_cmo_variables
    call ivo_cmo_read

    if (rank == 0) then
        call open_formatted_file(unit=unit_buf, file='BUF_write', status='replace', optional_action='write')
        write (unit_buf, '(6F22.16)') BUF(:)
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
                if (rank == 0) print *, 'all syminfo is zero, idx_irrep = ', idx_irrep
                nv0 = mo_end_idx - mo_start_idx + 1
                is_all_syminfo_zero = .true.
            else
                if (allocated(kappa)) then
                    nv0 = 0
                    do idx = mo_start_idx, mo_end_idx
                        ! If kappa is available, we check MJ instead of syminfo.
                        ! In atomic case, multiple D2h symmetries can have same MJ.
                        ! Ar(Atom) has linear symmetry, so we check Linear ID = abs(MJ) = 2*INDI - 1.
                        call atomic_id(syminfo(idx), kp, j, mj, ll)
                        indi = (abs(mj) + 1)/2
                        if (2*indi - 1 == isym_for_syminfo) then
                            nv0 = nv0 + 1
                        end if
                    end do
                else
                    nv0 = count(abs(syminfo(mo_start_idx:mo_end_idx)) == isym_for_syminfo)
                end if
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

    call ivo_cmo_write
    if (debug .and. rank == 0) print *, 'fockivo end'
    deallocate (BUF)
    deallocate (eval, syminfo)
contains
    subroutine create_itrfmo
        implicit none

        do iao = 1, num_ao
            do imo = 1, num_mo
                i0 = juck_up_idx + imo - 1
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
            if (is_all_syminfo_zero) then
                cnt = cnt + 1
                dmosym(cnt) = idx
            else
                if (allocated(kappa)) then
                    call atomic_id(syminfo(idx), kp, j, mj, ll)
                    indi = (abs(mj) + 1)/2
                    if (2*indi - 1 == isym_for_syminfo) then
                        cnt = cnt + 1
                        dmosym(cnt) = idx
                    end if
                else
                    if (abs(syminfo(idx)) == isym_for_syminfo) then
                        cnt = cnt + 1
                        dmosym(cnt) = idx
                    end if
                end if
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
end subroutine fockivo
